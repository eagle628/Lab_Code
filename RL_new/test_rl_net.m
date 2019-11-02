% general actor critic test
clear
close all

%% define model
% % % network version
seed1 = 8;
Node_number = 2;
net = network_swing_simple(Node_number, [1,2], [2,12]*1e-3, 1, [0.1,6], 0.8, seed1);
% net.Adj_ref = net.Adj_ref*0;
% local system
c_n = 1;
% set model
Ts = 0.1;
model = swing_network_model(net, c_n, Ts);
%% belief_N(What numbers do observe signal store ?)
belief_N = 3;
%% define init controller
seed_define = 2;
rng(seed_define)
% sys = ss(model.A, model.B, eye(2), [], model.Ts);
controller_n = 4;% state
controller_l = 1;% output
nnn = model.ny;
controller_m = belief_N*nnn;% input
A = diag(ones(1,(belief_N-1)*nnn),-nnn);
B = zeros(size(A,1),nnn);
B(1:nnn,1:nnn) = eye(nnn);
C = A;
D = B;
recorder = ss(A,B,C,D,model.Ts);
iter = 1;
disp('initial controller')
while true
    controller = drss(controller_n, controller_l, controller_m);
    loop = feedback(model.sys_local_discrete, controller*recorder, +1);
    if isstable(loop)
        break;
    end
    disp(iter)
    iter = iter+1;
end
% transform canon
recorder_sys = [A,B];

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

train_seed = 28;
Te = 50;
train.max_episode = 3000;
initial_set = zeros(model.nx, train.max_episode);
initial_set(1:end-model.rect_nx, :) = 2*rand(model.nx-model.rect_nx, train.max_episode)-1;
[x_all, rl_u_all, policy_snapshot, value_snapshot, reward_history] = train.train(initial_set, Te, train_seed);

%%
% test_ini = [-0.4;0];
% test_Te = 10;
% [x_all1, y_all1, rl_u_all1, t1, reward1] = train.sim(test_ini, test_Te);
% [x_all2, y_all2, mpc_u_all2, t2, reward2] = train.sim_lqrcontroller(test_ini, test_Te);
% 
%  figure
% plot(t1,y_all1)
% plot(t1,y_all1)
% grid on
% hold on
% plot(t1,y_all2,'--')
% norm(y_all1(1,:))
% norm(y_all2(1,:))
% 
% figure
% plot(t1, rl_u_all1)
% hold on
% grid on
% plot(t1, mpc_u_all2,'--')
