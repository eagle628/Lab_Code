% general actor critic test
clear
close all

% define model
l=0.5; M=0.15; g=9.8; eta=0.05; Ts=0.01;
model = pendulum_environment(l, M, g, eta, Ts);
% model = inverted_pendulum_model(l, M, g, eta, Ts);

%% Set value
basis_N = 11^2;
nnn = sqrt(basis_N);
range = [-4,4];
width = (range(2)-range(1))/(nnn-1);
m = (range(1):width:range(2))';
range = [-2,2];
width = (range(2)-range(1))/(nnn-1);
mm = (range(1):width:range(2))';
mu = [kron(m,ones(nnn,1)),repmat(mm,nnn,1)]; 
sigma = 0.25*ones(basis_N, 1);
RBF1 = Radial_Basis_Function(mu, sigma);
RBF2 = Radial_Basis_Function(mu, sigma);
sigma_pi = 0.5;
policy = Stocastic_Policy(RBF1, sigma_pi);
value  =  Value_base(RBF2);

opt_policy = TD_lambda(policy, 1e-4, 0.99);
opt_policy.target.pi_grad_enable = true;
opt_value = TD_lambda(value, 5e-4, 0.99);
%%
train = AC_episodic(model, opt_policy, opt_value);

train_seed = 28;
Te = 2;
train.max_episode = 3000;
initial_set = [rand(1, train.max_episode) - 0.5; zeros(1, train.max_episode)];
[x_all, rl_u_all, policy_snapshot, value_snapshot, reward_history] = train.train(initial_set, Te, train_seed);

%%
test_ini = [-0.4;0];
test_Te = 10;
[x_all1, y_all1, rl_u_all1, t1, reward1] = train.sim(test_ini, test_Te);
[x_all2, y_all2, mpc_u_all2, t2, reward2] = train.sim_lqrcontroller(test_ini, test_Te);

 figure
plot(t1,y_all1)
plot(t1,y_all1)
grid on
hold on
plot(t1,y_all2,'--')
norm(y_all1(1,:))
norm(y_all2(1,:))

figure
plot(t1, rl_u_all1)
hold on
grid on
plot(t1, mpc_u_all2,'--')
