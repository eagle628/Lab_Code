clear
close all

model = CRLMBC_test_model(0.5,0.15,9.8,0.05,0.01);

train = CRLMBC_test_train(model, 3, 121);

% figure
% [x, u_mpc, u_rl] = train.sim([0.4;0]);
% 
% figure
% plot(u_mpc);
% hold on
% plot(u_rl)

figure,
x = train.sim_noncontroller([0.4;0]);
plot(train.t,x)

figure,
x = train.sim_lqrcontroller([0.4;0]);
plot(train.t,x)
