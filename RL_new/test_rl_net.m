% general actor critic test
clear
close all

%% define model
% % % network version
seed1 = 1;
Node_number = 2;
net = network_swing_simple(Node_number, [1,2], [2,12]*1e-3, 1, [0.1,6], 0.8, seed1);
% net.Adj_ref = net.Adj_ref*0;
% local system
c_n = 1;
% set model
Ts = 0.1;
model = swing_network_model(net, c_n, Ts);
%% belief_N(What numbers do observe signal store ?)
belief_N = 1;
%% define init controller
seed_define = 3;
rng(seed_define)
% controller dimension
controller_n = 4;% state
controller_l = 1;% output
nnn = model.ny;
controller_m = belief_N*nnn;% input
% recorder sys define
A = diag(ones(1,(belief_N-1)*nnn),-nnn);
B = zeros(size(A,1),nnn);
B(1:nnn,1:nnn) = eye(nnn);
C = A;
D = B;
recorder = ss(A,B,C,D,model.Ts);
% explore
iter = 1;
disp('initial controller')
while true
    controller = drss(controller_n, controller_l, controller_m);
    loop = feedback(model.sys_local_discrete, controller*recorder, +1);
    if isstable(loop)
%         if order(minreal(loop,[],false)) == order(loop)
            break;
%         end
    end
    iter = iter+1;
end
fprintf('Found the intial controller for the %dth time.\n', iter)
% transform canon
recorder_sys = [A,B];
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
net_seed = 10;
py.deep_network_model.deep_model.reset_seed(int64(net_seed));
deep_net = py.deep_network_model.deep_model.simple_net();
% deep_net = py.deep_network_model.deep_model.simple_linear5_net();
gpu_id = int64(0);
if gpu_id >= 0
    deep_net.to_gpu(gpu_id);
end
% value = Chainer_Deep_Value(py.chainer.optimizers.SGD(pyargs('lr',0.01)).setup(deep_net), gpu_id);
value = Chainer_Deep_Value(py.chainer.optimizers.RMSprop(pyargs('lr',0.01)).setup(deep_net), gpu_id);
opt_value = Chainer_Deep_Optimizer(value);
% config
opt_value.constraint_enable = false;

%% set policy
ss_model = gen_ss_tridiag(controller_n,controller_m,controller_l);
ss_model.set_sys(controller);
apx_function2 = Dynamic_LTI_SS(ss_model);
sigma_pi = 10;
policy = Stocastic_Policy(apx_function2, sigma_pi);
opt_policy = TD_lambda(policy, 1e-9, 0);
opt_policy.constraint_enable = true;
opt_policy.target.pi_grad_enable = true;

%%
train = AC_episodic_for_net(model, opt_policy, opt_value, recorder_sys);

train_policy_seed = 28;
train_initial_seed = 1024;
Te = 50;
train.max_episode = 1000;
initial_set = zeros(model.nx, train.max_episode);
rng(train_initial_seed)
initial_set(1:end-model.rect_nx, :) = 2*rand(model.nx-model.rect_nx, train.max_episode)-1;
[x_all_train, u_all_train, policy_snapshot, value_snapshot, reward_train] = train.train(initial_set, Te, train_policy_seed);

%%
test_initial_seed = 256;
rng(test_initial_seed);
test_ini = zeros(model.nx, 1);
test_ini(2*model.c_n-1:2*model.c_n) = randn(model.local_nx ,1);
% test_ini(1:end-2) = randn(model.nx -2, 1);
test_Te = 1000;
[x_all_rl, y_all_rl, u_all_rl, t1, reward1_rl] = train.sim(test_ini, test_Te);
Q = eye(model.local_nx);
R = eye(model.nu);
[x_all_lqr, y_all_lqr, u_all_lqr, t2, reward_lqr] = train.sim_lqrcontroller(test_ini, test_Te, Q, R);
[x_all_elqr, y_all_elqr, u_all_elqr, t3, reward_elqr] = train.sim_extendlqrcontroller(test_ini, test_Te, controller_n-model.local_nx, Q, R);

figure('Name','Respose')
plot(t1, y_all_rl(2, :));
hold on, grid on
plot(t1, y_all_lqr(2, :), ':');
plot(t1, y_all_elqr(2, :), '--');

disp('RL')
disp(norm(y_all_rl(2, :)));
disp('LQR')
disp(norm(y_all_lqr(2, :)));
disp('Extend LQR')
disp(norm(y_all_elqr(2, :)));

%%
% savename = char(datetime);
% savename = strrep(savename,'/','-');
% savename = strrep(savename,':','-');
% savename = strrep(savename,' ','-');
% save(savename)
