clear
close all

rng(0)

model = CRLMBC_test_model(0.5,0.15,9.8,0.05,0.01);
model.set_observer_gain(0.5);

%% state feedback  and CRLMBC
% train = RL_state_feedback_and_observer_train(model, 3, 121);
train = RL_state_feedback_train(model, 3, 121);
% train = LQR_model_free_state_feedback_train(model, 10);
% K = train.train([0.4;0]);
% 
% y = train.sim([0.4;0], K);
% figure
% plot(train.t, y)

figure
[x, u_mpc, u_rl, omega] = train.actor_critic([0.4;0]);

figure
plot(u_mpc);
hold on
plot(u_rl)

%%
x_all = train.sim([0.4;0], omega);

figure
plot(train.t, x_all(:,1), 'r')
hold on

x = train.sim_lqrcontroller([0.4;0]);
plot(train.t,x(:,1),'b')

legend('RL+LQR','LQR')


%% IRL
% train = LQR_IRL_VI_alg_train(model, 10);
% train.train([0.4;0.2]);

% figure,
% plot(train.t, x);
