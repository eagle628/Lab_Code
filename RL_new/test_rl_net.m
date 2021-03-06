% general actor critic test
% clear
close all

%% define model
% % % network version
net_seed = 1;
Node_number = 30;
net = network_swing_simple(Node_number, [1,2], [2,12]*1e-3, 1, [0.1,6], 0.8, net_seed);
% net.Adj_ref = net.Adj_ref*0;
% local system
c_n = 1;
% set model
Ts = 0.1;
model = swing_network_model(net, c_n, Ts);
%% define belief sys
% belief_N(What numbers do observe signal store ?)
% belief_N = 20;
% recorder_sys = observation_accumulater(model.ny, belief_N);
% one step predictive pbserver
% model_rl_env = balred(model.RL_env_all({'yhat','what'},:), 6);
model_rl_env = model.sys_local_discrete;
recorder_sys = one_step_predictive_observer(model_rl_env.A, model_rl_env.B, model_rl_env.C);
%% define init controller
seed_define = 240;
rng(seed_define)
% controller dimension
controller_n = 4;% state
controller_l = 1;% output
controller_m = size(recorder_sys.C, 1);% input
% explore
iter = 1;
disp('initial controller')
while true
    controller = drss(controller_n, controller_l, controller_m);
    [a,b,c] = recorder_sys.connect(model.sys_local_discrete);
    loop = feedback(controller, ss(a,b,c,[],model.Ts), +1);
    if isstable(loop)
        break;
    end
    iter = iter+1;
end
fprintf('Found the intial controller for the %dth time.\n', iter)
% transform canon
figure('Name','Loop Pzmap')
pzmap(loop)


%% generate apx function for value
% % rbf etc.
% basis_N = 5;
% % generate line
% % basis function RBF
% range = [-5,5];
% width = (range(2)-range(1))/(basis_N-1);
% m = range(1):width:range(2);
% mu = m;
% for itr = 1 : belief_N*nnn - 1
% mu = combvec(mu, m); % col vector combinater function
% end
% mu = mu';
% sigma = 0.25*ones(size(mu, 1), 1);
% apx_function1 = Radial_Basis_Function(mu, sigma);
% value  =  Value_base(apx_function1);
% opt_value = TD_lambda(value, 5e-4, 0.99);
% deep net
deepnet_seed = 10;
py.deep_network_model.deep_model.reset_seed(int64(deepnet_seed));
% deep_net = py.deep_network_model.deep_model.simple_net();
% deep_net = py.deep_network_model.deep_model.simple_linear5_net();
% deep_net = py.deep_network_model.deep_model.linear3_net();
deep_net = py.deep_network_model.deep_model.linear5_net(pyargs('ratio',0));
gpu_id = int64(0);
if gpu_id >= 0
    deep_net.to_gpu(gpu_id);
end
% fixed target network enable
fixed_apx_function_enable = true;
% construct value
value_initial_lr = 1;
% chainer_optimizer = py.chainer.optimizers.SGD(pyargs('lr',initila_lr)).setup(deep_net);
chainer_optimizer = py.chainer.optimizers.RMSprop(pyargs('lr',value_initial_lr)).setup(deep_net);
value = Chainer_Deep_Value(chainer_optimizer, gpu_id, fixed_apx_function_enable);
opt_value = Chainer_Deep_Optimizer(value);
% config
opt_value.constraint_enable = false;
opt_value.trigger_enable = false;
opt_value.trigger_period = 1;
value_lr_inf_lim = 0.1;
opt_value.trigger_form = @(x) decay(x, value_lr_inf_lim, 1, (value_initial_lr-value_lr_inf_lim)/1000);

%% set policy
ss_model = gen_ss_tridiag(controller_n,controller_m,controller_l);
ss_model.set_sys(controller);
apx_function2 = Dynamic_LTI_SS(ss_model);
sigma_pi = 0.1;
policy = Stocastic_Policy(apx_function2, sigma_pi);
policy_initial_lr = 1e-1;
opt_policy = TD_lambda(policy, policy_initial_lr, 0);
% config
opt_policy.constraint_enable = true;
opt_policy.target.pi_grad_enable = false;
opt_policy.trigger_enable = false;
policy_lr_inf_lim = 1e-6;
opt_policy.trigger_period = 1;
opt_policy.trigger_form = @(x) decay(x, policy_lr_inf_lim, 1, (policy_initial_lr-policy_lr_inf_lim)/1000);

%% define train class
train = AC_episodic_for_net(model, opt_policy, opt_value, recorder_sys);
%% train contdition
train_policy_seed = 256;
train_initial_seed = 1024;
train_Te = 50;
train.max_episode = 2000;
train.snapshot = 100;
train.fixed_apx_function_period = 4;
train.gamma = 1;
train_initial_set = zeros(model.nx, train.max_episode);
rng(train_initial_seed)
% train_initial_set(1:end-model.rect_nx, :) = 2*rand(model.nx-model.rect_nx, train.max_episode)-1;
%% test condition
test_initial_seed = 64;
rng(test_initial_seed);
test_ini = zeros(model.nx, 1);
test_ini(2*model.c_n-1:2*model.c_n) = randn(model.local_nx ,1);
% test_ini(1:end-2) = randn(model.nx -2, 1);
test_Te = 100;
% initial 
[x_all_rl_initial, y_all_rl_initial, u_all_rl_initial, t, reward1_rl_initial] = train.sim(test_ini, test_Te);
z_all_rl_initial = train.model.evaluate(x_all_rl_initial);
%% Run train
figure('Name','Train Progress report')
[x_all_train, u_all_train, policy_snapshot, value_snapshot, history] = train.train(train_initial_set, train_Te, train_policy_seed);

%%

[x_all_rl_train, y_all_rl_train, u_all_rl_train, ~, reward1_rl_train] = train.sim(test_ini, test_Te);
z_all_rl_train = train.model.evaluate(x_all_rl_train);
[x_all_original, y_all_original, ~] = train.sim_original(test_ini, test_Te);
Q = eye(model.local_nx);
R = eye(model.nu);
[x_all_lqr, y_all_lqr, u_all_lqr, ~, reward_lqr] = train.sim_lqrcontroller(test_ini, test_Te, Q, R);
[x_all_elqr, y_all_elqr, u_all_elqr, ~, reward_elqr] = train.sim_extendlqrcontroller(test_ini, test_Te, controller_n-model.local_nx, Q, R);

figure('Name','Respose&Input')
subplot(1,2,1)
plot(t, z_all_rl_initial(2, :),'r:');
hold on, grid on
plot(t, z_all_rl_train(2, :),'m--');
plot(t, y_all_original(2, :),'k');
plot(t, y_all_lqr(2, :), 'b');
plot(t, y_all_elqr(2, :), 'g');
subplot(1,2,2)
plot(t, u_all_rl_initial,'r:')
hold on, grid on
plot(t, u_all_rl_train,'m--')
plot(t, u_all_lqr,'b');
plot(t, u_all_elqr,'g');

disp('Error-norm')
disp('RL-initial')
disp(norm(z_all_rl_initial(2, :)));
disp('RL-train')
disp(norm(z_all_rl_train(2, :)));
disp('Original')
disp(norm(y_all_original(2, :)));
disp('LQR')
disp(norm(y_all_lqr(2, :)));
disp('Extend-LQR')
disp(norm(y_all_elqr(2, :)));

disp('Reward')
disp('RL-initial')
disp(reward1_rl_initial);
disp('RL-train')
disp(reward1_rl_train);
disp('Original')
reward_original = cellfun(@(y)model.reward(y,0),...
                                    mat2cell(y_all_original, size(y_all_original,1), ones(1,size(y_all_original,2))));
tmp = 0:length(t)-1;
tmp = arrayfun(@(x)train.gamma^x,tmp);
reward_original = dot(reward_original, tmp);
disp(reward_original)
disp('LQR')
disp(reward_lqr);
disp('Extend-LQR')
disp(reward_elqr);

%%
savename = char(datetime);
savename = strrep(savename,'/','-');
savename = strrep(savename,':','-');
savename = strrep(savename,' ','-');
save(savename)

%% local
function out = decay(in, inf_lim, a, b)
    if in <= inf_lim
        out = in;
        return;
    end
    out = a*in-b;
end
