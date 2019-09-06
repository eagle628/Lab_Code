clear
close all

%% define model
% model = pendulum_model(0.5,0.15,9.8,0.05,0.01);
% % % network version
seed1 = 5;
Node_number = 3;
net = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,2], 0.1, seed1);
net.Adj_ref = net.Adj_ref*0;
% local system
c_n = 1;
% edge adjustment
seed2 = 10;
rng(seed2)
tmp_idx = find(net.Adj(c_n, :)~=0);
tmp_adj = 2 + 3*rand(1, length(tmp_idx));
net.Adj(c_n ,tmp_idx) = tmp_adj;
net.Adj(tmp_idx, c_n) = tmp_adj;
% net.plot()
% set model
Ts = 0.01;
model = swing_network_model(net, c_n, Ts);

%% belief_N
belief_N = 2;
%% define init controller
seed_define = 16;
rng(seed_define)
% sys = ss(model.A, model.B, eye(2), [], model.Ts);
controller_n = 2;% state
controller_l = 1;% output
controller_m = belief_N*model.ny;% input
A = diag(ones(1,belief_N-1), -1);
A = kron(A,repmat([1,0],model.ny,1));
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

% while true
%     controller = drss(2,1,2);
%     controller.Ts = model.Ts;
%     loop = loopsens(sys, controller);
%     if loop.Stable
%         break;
%     end
% end

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

%% generate apx function for policy
pi_sigma = 0.5;
apx_function = gen_ss_tridiag(controller_n,controller_m,controller_l);
apx_function.set_sys(controller);
policy = policy_dynamic(apx_function, pi_sigma);
% policy = policy_RBF(RBF1, pi_sigma);

%% trainer
Te = 10;
% train = general_actor_critic_with_eligibility_traces_episodic(model, policy, value, Te);
% train = general_actor_critic_with_eligibility_traces_episodic_on_POMDP(model, policy, value, Te, belief_N);
train = network_retro_by_actor_critic_with_eligibility_traces_episodic(model, policy, value, Te, belief_N);
train_seed = 28;
[x, u_mpc, u_rl, theta_mu_snapshot, theta_sigma_snapshot, w_snapshot, reward_history] = ...
    train.train([0.4, 0], train_seed);

%% plot
%%
test_ini = [0;0.4];
noise_seed = 10;
% x_all = train.sim(test_ini, [], 'parallel', mode_parallel);
x_rl = train.sim(test_ini,[],noise_seed);

figure
plot(train.t, x_rl(:,2), 'r')
hold on

x_lqr = train.sim_lqrcontroller(test_ini,[],noise_seed);
plot(train.t,x_lqr(:,2),'b')

x_ori = train.sim_original(test_ini,[],noise_seed);
plot(train.t,x_ori(:,2),'c')

legend('RL','LQR','ori')

% K = dlqr(model.A,model.B,train.Q,train.R);
% u = (-K*x')';
% 
% reward = 0;
% for itr = 2 : train.sim_N
%     reward = reward+ train.gamma^(itr-1)*train.reward(x(itr, :), u(itr-1, :));
%     
% end

%% test response
% % test_ini = [0;0];
% % noise_seed = 10;
% % x_lqr = train.sim_lqrcontroller(test_ini,[],noise_seed);
% % x_ori = train.sim_original(test_ini,[],noise_seed);
% % 
% % figure,plot(train.t,x_lqr(:,2),'r',train.t,x_ori(:,2),'b')
% % legend('lqr','original')
% % norm(x_lqr(:,2))
% % norm(x_ori(:,2))
