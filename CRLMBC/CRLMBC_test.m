clear
close all

rng(0)

model = CRLMBC_test_model(0.5,0.15,9.8,0.05,0.01);
% % % % % model.set_observer_gain(0.5);


%% state feedback  and CRLMBC
seed = 1024;
train = general_actor_critic_with_eligibility_traces_episodic(model, 3, 21^2, seed);

[x, u_mpc, u_rl, theta, w] = train.train([0.4, 0]);

figure
plot(u_mpc);
hold on
plot(u_rl)

%%
x_all = train.sim([0.4;0], theta);

figure
plot(train.t, x_all(:,1), 'r')
hold on

x = train.sim_lqrcontroller([0.4;0]);
plot(train.t,x(:,1),'b')

legend('RL+LQR','LQR')


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
