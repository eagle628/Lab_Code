clear
close all

%% define model
model = pendulum_model(0.5,0.15,9.8,0.05,0.01);

%% define init controller
sys = ss(model.A, model.B, eye(2), [], model.Ts);
controller_n = 2;% state
controller_l = 1;% output
controller_m = 2;% input
while true
    controller = drss(2,1,2);
    controller.Ts = model.Ts;
    loop = loopsens(sys, controller);
    if loop.Stable
        break;
    end
end

%% generate apx function for value
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
value  =  value_RBF(RBF1);

%% generate apx function for policy
apx_function = gen_ss_tridiag(controller_n,controller_m,controller_l);
apx_function.set_sys(controller);
pi_sigma = 0.5;
policy = policy_dynamic(apx_function, pi_sigma);

%% trainer
Te = 2;
train = general_actor_critic_with_eligibility_traces_episodic(model, policy, value, Te);
mode_parallel = 'off';
train_seed = 28;
[x, u_mpc, u_rl, theta_mu_snapshot, theta_sigma_snapshot, w_snapshot, reward_history] = ...
    train.train([0.4, 0], train_seed, 'parallel', mode_parallel);

%% plot
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
