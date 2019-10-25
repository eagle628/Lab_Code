clear
close all

%% define model
% % % network version
seed1 = 8;
Node_number = 2;
net = network_swing_simple(Node_number, [1,2], [2,12]*1e-3, 1, [0.1,6], 0.8, seed1);
net.Adj_ref = net.Adj_ref*0;
% local system
c_n = 1;
% % edge adjustment
% seed2 = 10;
% rng(seed2)
% tmp_idx = find(net.Adj(c_n, :)~=0);
% tmp_adj = 5 + 1*rand(1, length(tmp_idx));
% net.Adj(c_n ,tmp_idx) = tmp_adj;
% net.Adj(tmp_idx, c_n) = tmp_adj;
% net.plot()
% set model
Ts = 0.01;
model = swing_network_model(net, c_n, Ts);
disp('IsStable System')
disp(isstable(c2d(model.sys_all, Ts)))
%% belief_N(What numbers do observe signal store ?)
belief_N = 2;
%% define init controller
seed_define = 2;
rng(seed_define)
% sys = ss(model.A, model.B, eye(2), [], model.Ts);
controller_n = 4;% state
controller_l = 1;% output
nnn = model.ny+model.nw;
controller_m = belief_N*nnn;% input
A = diag(ones(1,(belief_N-1)*nnn),-nnn);
B = zeros(size(A,1),nnn);
B(1:nnn,1:nnn) = eye(nnn);
C = A;
D = B;
recorder = ss(A,B,C,D,model.Ts);
P = model.sys_local({'y','w'},{'u'});
% P.C = zeros(size(P.C));
% P.C(1:model.ny,1:model.ny) = eye(model.ny);
P = c2d(P, model.Ts);
iter = 1;
disp('initial controller')
while true
    controller = drss(controller_n, controller_l, controller_m);
    loop = loopsens(P, -controller*recorder);
    if loop.Stable
        break;
    end
    disp(iter)
    iter = iter+1;
end
% transform canon
recorder_sys = [A,B];

%% generate apx function for value
% % parmeter
% basis_N = 11;
% % generate line
% % basis function RBF
% range = [-1,1];
% width = (range(2)-range(1))/(basis_N-1);
% m = range(1):width:range(2);
% mu = m;
% for itr = 1 : belief_N*nnn - 1
% mu = combvec(mu, m); % col vector combinater function
% end
% mu = mu';
% sigma = 0.25*ones(size(mu, 1), 1);
% RBF1 = Radial_Basis_Function(size(mu, 1), mu, sigma);
% value  =  value_RBF(RBF1);
% value_update_rule = TD_lambda(1e-5, 0.999);
% opt_value = optimizer(value_update_rule, value);
% % deep ver
net_seed = 10;
py.deep_network_model.deep_model.reset_seed(int64(net_seed));
deep_net = py.deep_network_model.deep_model.simple_net();
gpu_id = int64(0);
if gpu_id >= 0
    deep_net.to_gpu(gpu_id);
end
value = value_chainer_deep_net(py.chainer.optimizers.SGD(pyargs('lr',0.01)).setup(deep_net));
value_update_rule = deep_update_rule();
opt_value = optimizer(value_update_rule, value);

%% generate apx function for policy
pi_sigma = 1;
% apx_function = gen_tf(controller_n,controller_m,controller_l);
apx_function = gen_ss_tridiag(controller_n,controller_m,controller_l);
apx_function.set_sys(controller);
% policy = policy_dynamic_tf(apx_function, pi_sigma);
policy = policy_dynamic_ss(apx_function, pi_sigma);
policy_update_rule = TD_lambda(1e-5, 0);
opt_policy = optimizer(policy_update_rule, policy);

%% trainer
train = network_retro_by_AC_episodic(model, opt_policy, opt_value, recorder_sys);
train.max_episode = 10;
train_seed = 28;

Te = 5;
mode_parallel = false;
[x, u_mpc, u_rl, policy_snapshot, value_snapshot, reward_history] = ...
    train.train([0.4, 0], Te, train_seed,'mode-parrallerl', mode_parallel);

% load test_rl_net30_1010.mat obj
% train.opt_policy.approximate_function_class.set_params(obj.opt_policy.approximate_function_class.get_params())
% train.opt_policy.approximate_function_class.set_params(obj.opt_value.approximate_function_class.get_params())
%% plot
test_ini = [0;0.4];
test_Te  = 10;
t = (0:model.Ts:test_Te)';
noise_seed = 1;
x_rl = train.sim(test_ini, test_Te, [], [], noise_seed,'mode-parrallerl', mode_parallel);

figure
plot(t, x_rl(:,2), 'r')
hold on

x_lqr = train.sim_lqrcontroller(test_ini,test_Te,[],noise_seed);
plot(t,x_lqr(:,2),'b')

x_ori = train.sim_original(test_ini,test_Te,[],noise_seed);
plot(t,x_ori(:,2),'c')

net.add_controller(model.c_n, balred(model.sys_env, controller_n-2), train.Q_c, train.R)
con_sys = net.get_sys_controlled(model.sys_all);
rng(noise_seed)
ddd = randn(length(t), 2);
con_sys = con_sys({strcat('y_node',num2str(model.c_n))},{strcat('d_node',num2str(model.c_n))});
x0 = zeros(order(con_sys), 1);
x0(1:2) = test_ini;
x_compared = lsim(con_sys, ddd, t, x0);
net.controllers = {};
plot(t,x_compared(:,2),'k:')
rng('shuffle')

legend('RL','LQR','ori','compared')

disp('RL')
disp(norm(x_rl(:,2)))
disp('LQR')
disp(norm(x_lqr(:,2)))
disp('Original')
disp(norm(x_ori(:,2)))
disp('compared')
disp(norm(x_compared(:,2)))

%% parfor test
% % % % test_ini = [0;0.4];
% % % % test_Te  = 10;
% % % % t = (0:model.Ts:test_Te)';
% % % % parfor_N = 1e3;
% % % % 
% % % % % extend retro
% % % % net.add_controller(model.c_n, balred(model.sys_env, controller_n-2), train.Q_c, train.R)
% % % % con_sys = net.get_sys_controlled(model.sys_all);
% % % % con_sys = con_sys({strcat('y_node',num2str(model.c_n))},{strcat('d_node',num2str(model.c_n))});
% % % % net.controllers = {};
% % % % x0 = zeros(order(con_sys), 1);
% % % % x0(1:2) = test_ini;
% % % % % get memory
% % % % norm_rl_set = zeros(parfor_N, 1);
% % % % norm_lqr_set = zeros(parfor_N, 1);
% % % % norm_ori_set = zeros(parfor_N, 1);
% % % % norm_compared_set = zeros(parfor_N, 1);
% % % % % loop
% % % % parfor_progress(parfor_N)
% % % % parfor noise_seed = 1 : parfor_N
% % % %     x_rl = train.sim(test_ini, test_Te, [], [], noise_seed,'mode-parrallerl', mode_parallel);
% % % %     norm_rl_set(noise_seed, 1) = norm(x_rl(:,2));
% % % %     
% % % %     x_lqr = train.sim_lqrcontroller(test_ini,test_Te,[],noise_seed);
% % % %     norm_lqr_set(noise_seed, 1) = norm(x_lqr(:,2));
% % % %     
% % % %     x_ori = train.sim_original(test_ini,test_Te,[],noise_seed);
% % % %     norm_ori_set(noise_seed, 1) = norm(x_ori(:,2));
% % % %     
% % % %     rng(noise_seed)
% % % %     ddd = randn(length(t), 2);
% % % %     x_compared = lsim(con_sys, ddd, t, x0);
% % % %     norm_compared_set(noise_seed, 1) = norm(x_compared(:,2));
% % % %     rng('shuffle')
% % % %     
% % % %     parfor_progress();
% % % % end
% % % % parfor_progress(0);
% % % % 
% % % % figure
% % % % boxplot([norm_ori_set,norm_lqr_set,norm_rl_set,norm_compared_set],'Labels',{'ori','lqr','rl','compared'})
%%
% mail_message('End')

% savename = char(datetime);
% savename = strrep(savename,'/','-');
% savename = strrep(savename,':','-');
% savename = strrep(savename,' ','-');
% save(savename)
%% local
