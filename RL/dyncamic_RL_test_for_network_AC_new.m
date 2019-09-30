clear
close all

%% define model
% % % network version
seed1 = 8;
Node_number = 2;
net = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,2], 0.8, seed1);
net.Adj_ref = net.Adj_ref*0;
% local system
c_n = 1;
% edge adjustment
seed2 = 10;
rng(seed2)
tmp_idx = find(net.Adj(c_n, :)~=0);
tmp_adj = 5 + 1*rand(1, length(tmp_idx));
net.Adj(c_n ,tmp_idx) = tmp_adj;
net.Adj(tmp_idx, c_n) = tmp_adj;
% net.plot()
% set model
Ts = 0.01;
model = swing_network_model(net, c_n, Ts);

%% belief_N(What numbers do observe signal store ?)
belief_N = 2;
%% define init controller
seed_define = 2;
rng(seed_define)
% sys = ss(model.A, model.B, eye(2), [], model.Ts);
controller_n = 2;% state
controller_l = 1;% output
controller_m = belief_N*model.ny;% input
A = diag(ones(1,(belief_N-1)*model.ny),-model.ny);
B = zeros(size(A,1),model.ny);
B(1:model.ny,1:model.ny) = eye(model.ny);
C = A;
D = B;
recorder = ss(A,B,C,D,model.Ts);
P = model.sys_local({'y'},{'u'});
P = c2d(P, model.Ts);
iter = 1;
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
% parmeter
basis_N = 11;
% generate line
% basis function RBF
range = [-1,1];
width = (range(2)-range(1))/(basis_N-1);
m = range(1):width:range(2);
mu = m;
for itr = 1 : belief_N*model.ny - 1
mu = combvec(mu, m); % col vector combinater function
end
mu = mu';
sigma = 0.25*ones(size(mu, 1), 1);
RBF1 = Radial_Basis_Function(size(mu, 1), mu, sigma);
value  =  value_RBF(RBF1);
value_update_rule = TD_lambda();
opt_value = optimizer(value_update_rule, value);

%% generate apx function for policy
pi_sigma = 1;
apx_function = gen_tf(controller_n,controller_m,controller_l);
apx_function.set_sys(controller);
policy = policy_dynamic_tf(apx_function, pi_sigma);
policy_update_rule = TD_lambda();
opt_policy = optimizer(policy_update_rule, policy);

%% trainer
train = network_retro_by_AC_episodic(model, opt_policy, opt_value, recorder_sys);
train_seed = 28;

Te = 5;
mode_parallel = true;
[x, u_mpc, u_rl, policy_snapshot, value_snapshot, reward_history] = ...
    train.train([0.4, 0], Te, train_seed,'mode-parrallerl', mode_parallel);

% train.train([0.4, 0], train_seed, 'Invalid-Constraint', true);

%% plot
%%
test_ini = [0;0.4];
test_Te  = 10;
t = (0:model.Ts:test_Te)';
noise_seed = 100;
 
x_rl = train.sim(test_ini, test_Te, [], [], noise_seed,'mode-parrallerl', mode_parallel);

figure
plot(t, x_rl(:,2), 'r')
hold on

x_lqr = train.sim_lqrcontroller(test_ini,test_Te,[],noise_seed);
plot(t,x_lqr(:,2),'b')

x_ori = train.sim_original(test_ini,test_Te,[],noise_seed);
plot(t,x_ori(:,2),'c')

legend('RL','LQR','ori')

disp('RL')
disp(norm(x_rl(:,2)))
disp('LQR')
disp(norm(x_lqr(:,2)))
disp('Original')
disp(norm(x_ori(:,2)))

%% local
