% function test_rl_net_function(data)
clear
close all

%%
data.net_seed = 8;
data.node_number = 6;
data.net_local = 1;
data.Ts = 0.1;
data.belief_N = 32;
data.initial_controller_seed = 244;
data.controller_state_number = 4;
data.value_deep_net_seed = 10;
data.value_fixed_apx_funciton_enable = false;
data.value_initial_lr = 0.01;
% data.value_initial_lr = 2e-5;
data.value_constraint_enable = false;
data.value_trigger_enable = false;
data.value_trigger_period = 1;
data.value_lr_inf_lim = 0.1;
data.policy_pi_sigma = 1000;
data.policy_initial_lr = 1e-7; %% for AC
% data.policy_initial_lr = 1e-12; %% for REINFORCE
data.policy_constraint_enable = true;
data.policy_pi_grad_enable = false;
data.policy_trigger_enable = false;
data.policy_lr_inf_lim = 1e-2;
data.policy_trigger_period = 1;
data.train_seed = 28;
data.train_initial_seed = 1024;
data.train_Te = 50;
data.train_max_episode = 12000;
data.train_snapshot = 100;
data.train_fixed_apx_fucntion_period = 1;
data.train_gamma = 1;
data.train_render_enable = false;
data.test_initial_seed = 64;
data.test_Te = 50;

%% define model
% % % network version
net = network_swing_simple(data.node_number, [1,2], [2,12]*1e-3, 1, [0.1,6], 0.8, data.net_seed);
net.Adj_ref = net.Adj_ref*0;
% set model
% [~, sys_env] = net.get_sys_local(data.net_local);
% model = swing_network_model(net, data.net_local, data.Ts, c2d(sys_env,data.Ts));
model = swing_network_model(net, data.net_local, data.Ts);
%% define belief sys
% belief_N(What numbers do observe signal store ?)
recorder_sys = observation_accumulater(model.ny, data.belief_N);
% one step predictive pbserver
% model_rl_env = balred(model.RL_env_all({'yhat','what'},:), 6);
% recorder_sys = one_step_predictive_observer(model_rl_env.A, model_rl_env.B, model_rl_env.C);
% recorder_sys = one_step_predictive_observer(model.dynamics_A, zeros(size(model.dynamics_B)), model.dynamics_C);
%% define init controller
% rng(data.initial_controller_seed)
% iter = 0;
% while true
%     controller = randn(size(recorder_sys.D, 1),1);
%     controller = Static_Gain(size(recorder_sys.D, 1), model.nu, controller);
%     [a2,b2,c2,d2] = controller.get_ss();
%     [a,b,c] = recorder_sys.connect(model.sys_local_discrete);
%     loop = feedback(ss(a2,b2,c2,d2,model.Ts), ss(a,b,c,[],model.Ts), +1);
%     if isstable(loop)
%         break;
%     end
%     iter = iter+1;
% end
% fprintf('Found the intial controller for the %dth time.\n', iter)

% init dynamic controller
rng(data.initial_controller_seed)
% controller dimension
controller_l = model.nu;% output
controller_m = size(recorder_sys.C, 1);% input
% explore
iter = 1;
disp('initial controller')
while true
    controller = drss(data.controller_state_number, controller_l, model.ny);
    [a,b,c,d] = ssdata(controller);
    % % % 
% % % % %     ss_model = gen_ss_tridiag(data.controller_state_number-2, 1, 1);
% % % % %     ss_model2 = gen_ss_Eretro_controller(ss_model, c2d(model.sys_local,model.Ts,model.discrete_type));
% % % % %     ss_model2.set_params(randn(1,ss_model2.N))
% % % % %     [a,b,c,d] = ss_model2.get_ss();
    controller = ss(a,...
                    [b,zeros(size(b,1),size(recorder_sys.C, 1)-size(b,2))],...
                    c,...
                    [d,zeros(size(d,1),size(recorder_sys.C,1)-size(d,2))],model.Ts);
    [a,b,c] = recorder_sys.connect(model.sys_local_discrete);
    loop = feedback(controller, ss(a,b,c,[],model.Ts), +1);
    if isstable(loop)
        break;
    end
    iter = iter+1;
end
fprintf('Found the intial controller for the %dth time.\n', iter)
% transform canon
% figure('Name','Loop Pzmap')
% pzmap(loop)
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
py.deep_network_model.deep_model.reset_seed(int64(data.value_deep_net_seed));
% deep_net = py.deep_network_model.deep_model.simple_net();
% deep_net = py.deep_network_model.deep_model.simple_linear5_net();
deep_net = py.deep_network_model.deep_model.linear3_net(pyargs('n_mid_units',int64(10)));
% deep_net = py.deep_network_model.deep_model.linear5_net(pyargs('ratio',0));
gpu_id = int64(0);
if gpu_id >= 0
    deep_net.to_gpu(gpu_id);
end
% construct value
% chainer_optimizer = py.chainer.optimizers.Adam(pyargs('alpha',data.value_initial_lr)).setup(deep_net);
chainer_optimizer = py.chainer.optimizers.SGD(pyargs('lr',data.value_initial_lr)).setup(deep_net);
% chainer_optimizer = py.chainer.optimizers.RMSprop(pyargs('lr',data.value_initial_lr)).setup(deep_net);
value = Chainer_Deep_Value(chainer_optimizer, gpu_id, data.value_fixed_apx_funciton_enable);
opt_value = Chainer_Deep_Optimizer(value);
% config
opt_value.constraint_enable = data.value_constraint_enable;
opt_value.trigger_enable = data.value_trigger_enable;
opt_value.trigger_period = data.value_trigger_period;
value_lr_inf_lim = data.value_lr_inf_lim;
opt_value.trigger_form = @(x) decay(x, value_lr_inf_lim, 1, (value_initial_lr-value_lr_inf_lim)/1000);
%% set policy
ss_model = gen_ss_tridiag(data.controller_state_number, controller_m, controller_l);
ss_model.set_sys(controller);
initial_policy_params = ss_model.get_params();
% % % % % Store = load('C:\Users\NaoyaInoue\Documents\GitHub\Lab_Code_new\data_test_global\controller_disable\many_dataset\2019-12-10=15-05-29.mat');
% % % % % [a,b,c,d] = Store.ss_model.get_ss(Store.x_pso);
% % % % % controller = ss(a,...
% % % % %                 [b,zeros(size(b,1),size(recorder_sys.C, 1)-size(b,2))],...
% % % % %                 c,...
% % % % %                 [d,zeros(size(d,1),size(recorder_sys.C,1)-size(d,2))],model.Ts);
% % % % % ss_model.set_sys(controller);
% load tmp x_test
% ss_model.set_params(x_test');
% % % % ss_model = gen_ss_Eretro_controller(ss_model, c2d(model.sys_local,model.Ts,model.discrete_type));
% % % % ss_model.set_params(ss_model2.get_params());
apx_function2 = Dynamic_LTI_SS(ss_model);
% apx_function2 = Static_Gain(size(recorder_sys.D, 1), model.nu, controller.theta);
policy = Stocastic_Policy(apx_function2, data.policy_pi_sigma);
policy_initial_lr = data.policy_initial_lr;
opt_policy = TD_lambda(policy, policy_initial_lr, 0, 1);
% config
opt_policy.constraint_enable = data.policy_constraint_enable;
opt_policy.target.pi_grad_enable = data.policy_pi_grad_enable;
opt_policy.trigger_enable = data.policy_trigger_enable;
policy_lr_inf_lim = data.policy_lr_inf_lim;
opt_policy.trigger_period = data.policy_trigger_period;
opt_policy.trigger_form = @(x) decay(x, policy_lr_inf_lim, 1, (policy_initial_lr-policy_lr_inf_lim)/1000);
%% define train class
train = AC_episodic_for_net(model, opt_policy, opt_value, recorder_sys);
% train = REINFORCE_MC_episodic_for_net(model, opt_policy, recorder_sys);
%% train contdition
train.render_enable = data.train_render_enable;
train.max_episode = data.train_max_episode;
train.snapshot = data.train_snapshot;
train.fixed_apx_function_period = data.train_fixed_apx_fucntion_period; %% For AC
train.value_pretraining_period = 1000; % For AC
train.gamma = data.train_gamma;
train_initial_set = zeros(model.nx, train.max_episode);
rng(data.train_initial_seed)
% train_initial_set(1:end-model.rect_nx, :) = 2*rand(model.nx-model.rect_nx, train.max_episode)-1;
%% Run train
if data.train_render_enable
    figure('Name','Train Progress report')
end
[x_all_train, u_all_train, policy_snapshot, value_snapshot, history] = train.train(train_initial_set, data.train_Te, data.train_seed);
%% test
test_N = 100;
reward_original_set = zeros(test_N, 1);
reward_lqr_set = zeros(test_N, 1);
reward_elqr_set = zeros(test_N, 1);
reward_rl_initial_set = zeros(test_N, 1);
reward_rl_train_set = zeros(test_N, 1);
rng(data.test_initial_seed);
test_ini = zeros(model.nx, 1);
% test_ini(2*model.c_n-1:model.rect_nx) = randn(model.rect_nx ,1);
test_Te = data.test_Te;
for test_d_seed = 1:test_N
% test_d_seed = 28;
% initial
final_policy_params = train.opt_policy.target.apx_function.form.get_params();
train.opt_policy.target.apx_function.form.set_params(initial_policy_params);
[x_all_rl_initial, y_all_rl_initial, u_all_rl_initial, t, reward1_rl_initial] = train.sim(test_ini, test_Te, test_d_seed);
z_all_rl_initial = train.model.evaluate(x_all_rl_initial);

% final
train.opt_policy.target.apx_function.form.set_params(final_policy_params);
[x_all_rl_train, y_all_rl_train, u_all_rl_train, ~, reward1_rl_train] = train.sim(test_ini, test_Te, test_d_seed);
z_all_rl_train = train.model.evaluate(x_all_rl_train);
[x_all_original, y_all_original, ~] = train.sim_original(test_ini, test_Te, test_d_seed);reward_original = cellfun(@(y)model.reward(y,0),...
                                    mat2cell(y_all_original, size(y_all_original,1), ones(1,size(y_all_original,2))));
tmp = 0:length(t)-1;
tmp = arrayfun(@(x)train.gamma^x,tmp);
reward_original = dot(reward_original, tmp);
Q = diag([1,1]);
R = eye(model.nu);
test_ini_lqr = test_ini(1:order(model.sys_all)+model.local_nx);
[x_all_lqr, y_all_lqr, u_all_lqr, ~, reward_lqr] = train.sim_lqrcontroller(test_ini_lqr, test_Te, test_d_seed, Q, R);
[x_all_elqr, y_all_elqr, u_all_elqr, ~, reward_elqr] = train.sim_extendlqrcontroller(test_ini_lqr, test_Te, data.controller_state_number-model.local_nx, test_d_seed, Q, R);


reward_original_set(test_d_seed) = reward_original;
reward_lqr_set(test_d_seed) = reward_lqr;
reward_elqr_set(test_d_seed) = reward_elqr;
reward_rl_initial_set(test_d_seed) = reward1_rl_initial;
reward_rl_train_set(test_d_seed) = reward1_rl_train;

end

figure
boxplot(...
    [reward_rl_initial_set,reward_rl_train_set,reward_lqr_set,reward_elqr_set],...
    {'Initial','Train','SLQR','ELQR'}...
    )
ylim([-inf 0])
f = gcf;
f.Children.YGrid = 'on';
f.Children.YMinorGrid = 'on';
f.Children.FontName = 'SanSelif';
f.Children.FontSize = 16;
f.Units = 'centimeters';
f.Position = [0,0,15,10.5];
xlabel('制御器設計方法')
ylabel('収益')
title('制御器設計手法による収益の分布')

% figure('Name','Respose&Input')
% subplot(3,1,1)
% plot(t, z_all_rl_initial(1, :),'r:');
% hold on, grid on
% plot(t, z_all_rl_train(1, :),'m--');
% plot(t, y_all_original(1, :),'k');
% plot(t, y_all_lqr(1, :), 'b');
% plot(t, y_all_elqr(1, :), 'g');
% subplot(3,1,2)
% plot(t, z_all_rl_initial(2, :),'r:');
% hold on, grid on
% plot(t, z_all_rl_train(2, :),'m--');
% plot(t, y_all_original(2, :),'k');
% plot(t, y_all_lqr(2, :), 'b');
% plot(t, y_all_elqr(2, :), 'g');
% subplot(3,1,3)
% plot(t, u_all_rl_initial,'r:')
% hold on, grid on
% plot(t, u_all_rl_train,'m--')
% plot(t, u_all_lqr,'b');
% plot(t, u_all_elqr,'g');
% 
% disp('Error-norm')
% disp('RL-initial')
% disp(norm(z_all_rl_initial(2, :)));
% disp('RL-train')
% disp(norm(z_all_rl_train(2, :)));
% disp('Original')
% disp(norm(y_all_original(2, :)));
% disp('LQR')
% disp(norm(y_all_lqr(2, :)));
% disp('Extend-LQR')
% disp(norm(y_all_elqr(2, :)));
% 
% disp('Reward')
% disp('RL-initial')
% disp(reward1_rl_initial);
% disp('RL-train')
% disp(reward1_rl_train);
% disp('Original')
% disp(reward_original)
% disp('LQR')
% disp(reward_lqr);
% disp('Extend-LQR')
% disp(reward_elqr);
%%
% savename = DataStruct2FileName(data);
savetime = char(datetime);
savetime = strrep(savetime,'/','-');
savetime = strrep(savetime,':','-');
savetime = strrep(savetime,' ','=');
save(savetime)
%% local
function out = decay(in, inf_lim, a, b)
    if in <= inf_lim
        out = in;
        return;
    end
    out = a*in-b;
end

% end
