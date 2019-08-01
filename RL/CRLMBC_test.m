% general actor critic test
clear
close all

rng(0)

model = CRLMBC_test_model(0.5,0.15,9.8,0.05,0.01);
% % % % % model.set_observer_gain(0.5);


%% state feedback  and CRLMBC
Te = 2;
% generate apx function
basis_N = 11^2;
nnn = sqrt(basis_N);
range = [-2,2];
width = (range(2)-range(1))/(nnn-1);
m = (range(1):width:range(2))';
range = [-4,4];
width = (range(2)-range(1))/(nnn-1);
mm = (range(1):width:range(2))';
mu = [kron(m,ones(nnn,1)),repmat(mm,nnn,1)]; 
sigma = 0.25*ones(basis_N, 1);
RBF1 = Radial_Basis_Function(basis_N, mu, sigma);
sigma_pi = sqrt(1);
policy = policy_RBF(RBF1, sigma_pi);
value  =  value_RBF(RBF1);
% set train
train = general_actor_critic_with_eligibility_traces_episodic(model, policy, value, Te);

mode_parallel = 'off';
train_seed = 28;
[x, u_mpc, u_rl, theta_mu_snapshot, theta_sigma_snapshot, w_snapshot, reward_history, F] = train.train([0.4, 0], train_seed, 'parallel', mode_parallel);

figure
plot(u_mpc);
hold on
plot(u_rl)

%%
x_all = train.sim([0.4;0], [], 'parallel', mode_parallel);

figure
plot(train.t, x_all(:,1), 'r')
hold on

x = train.sim_lqrcontroller([0.4;0]);
plot(train.t,x(:,1),'b')

legend('RL','LQR')


K = dlqr(model.A,model.B,train.Q,train.R);
u = (-K*x')';

reward = 0;
for itr = 2 : train.sim_N
    reward = reward+ train.gamma^(itr-1)*train.reward(x(itr, :), u(itr-1, :));
end

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
