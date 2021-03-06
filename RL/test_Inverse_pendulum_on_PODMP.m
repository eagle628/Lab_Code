% general actor critic test
clear
close all

rng(0)

% model = inverted_pendulummodel(0.5,0.15,9.8,0.05,0.01);
% model = inverted_pendulummodel(0.5,0.15,9.8,0,0.01);
model = pendulum_model(0.5,0.15,9.8,0.05,0.01);
% % % % % model.set_observer_gain(0.5);


%% state feedback  and CRLMBC
Te = 2;
% generate apx function
% parmeter
belief_N = 2;
basis_N = 11;
% generate line
% basis function RBF
range = [-1,1];
width = (range(2)-range(1))/(basis_N-1);
m = range(1):width:range(2);
mu = m;
for itr = 1 : belief_N * model.ny-1
mu = combvec(mu, m); % col vector combinater function
end
mu = mu';
sigma = 0.25*ones(size(mu, 1), 1);
RBF1 = Radial_Basis_Function(size(mu, 1), mu, sigma);
% basis funciton polynomial
% RBF1 = polinomial_basis_function(basis_N^(belief_N * model.ny), mu);
sigma_pi = 0.5;
policy = policy_RBF(RBF1, sigma_pi);
value  =  value_RBF(RBF1);
% set train
train = general_actor_critic_with_eligibility_traces_episodic_on_POMDP(model, policy, value, Te, belief_N);

Input_Clipping = 10;
TD_Error_Clipping = 5;
Fix_Target_Network = 100;
train_seed = 28;

% S = load('test.mat');
% train.policy.set_params(S.theta_mu_snapshot(:, end));
% train.policy.set_policy_sigma(S.theta_sigma_snapshot(:, end));
% train.value.set_params(S.w_snapshot(:, end));

[x, u_rl, theta_mu_snapshot, theta_sigma_snapshot, w_snapshot, reward_history, F1] = ...
    train.train([0.4, 0], train_seed, 'Input-Clipping',Input_Clipping);


% ,'TD-Error-Clipping',TD_Error_Clipping
% ,'Fix-Target-Network',Fix_Target_Network

% figure
% plot(u_mpc);
% hold on
% plot(u_rl)

%%
x_all = train.sim([-0.01;0], []);

figure
plot(train.t, x_all(:,1), 'r')
hold on

% x = train.sim_lqrcontroller([0.4;0]);
% plot(train.t,x(:,1),'b')
% 
% legend('RL','LQR')

% 
% K = dlqr(model.A,model.B,train.Q,train.R);
% u = (-K*x')';
% 
% reward = 0;
% for itr = 2 : train.sim_N
%     reward = reward+ train.gamma^(itr-1)*train.reward(x(itr, :), u(itr-1, :));
% end

%% apx
% % % % rng(0)
% % % % 
% % % % model = test_apx_model();
% % % % model.set_observer_gain(1);
% % % % 
% % % % train = RL_state_feedback_and_observer_train(model, 3, 121);
% % % % 
% % % % figure
% % % % [true_x, apx_x, u_mpc, u_rl, omega] = train.actor_critic(ones(model.true_nx,1));
