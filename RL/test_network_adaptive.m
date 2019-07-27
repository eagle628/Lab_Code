close all
clear

%% Release GPU memory
gpuDevice(1);

%%

seed1 = 6;
Node_number = 4;
net = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,2], 0.1, seed1);
net.Adj_ref = net.Adj_ref*0;

% local system
c_n = 1;
% edge adjustment
tmp_idx = find(net.Adj(c_n, :)~=0);
tmp_adj = 2 + 3*rand(1, length(tmp_idx));
net.Adj(c_n ,tmp_idx) = tmp_adj;
net.Adj(tmp_idx, c_n) = tmp_adj;

% train condition
Ts = 0.01;
Te = 1;

% model = swing_network_model(net, c_n, Ts);
[sys_local, sys_env] = net.get_sys_local(c_n);
apx_nx = 6;%net.N*2 -2;
% netwrok identificaiton
sys_all = net.get_sys();
for nnn = 1 : net.N
    sys_all = net.add_io(sys_all, nnn, strcat('node',num2str(nnn)));
end


%% Sngle Test
model = swing_network_model(net, c_n, Ts);

train_seed = 1024;
belief_N = 40;
train = netwrok_retro_by_adaptive_and_POMDP_episodic(model, belief_N, Te);

initial = [0,0.1];
%     initial = [0;1;];
[local_x_all, rl_u_all, theta_snpashot, reward_history] = train.train(initial, train_seed);
    %%
% %     sim_initial = [0;0;];
% %     [local_x_all1, env_x_all1, rect_x_all1, y_xhat_w_v_all1, rl_u_all1] = train.sim([], sim_initial, 6);
% %     [local_x_all2, env_x_all2, rect_x_all2, y_xhat_w_v_all2] = train.sim_lqrcontroller(sim_initial, 6);
    
%     save(strcat('RL_net',num2str(Node_number),'_discount99_', name{idx}, '_reward_fix_Q_episode_1000_No2'));
%     save(strcat('RL_net',num2str(Node_number),'_discount99_', name{idx}, '_reward_fix_Q_episode_2000_NonLQR'));
%     save(strcat('RL_net',num2str(Node_number),'_discount99_', name{idx}, '_reward_only_local_episode_',num2str(train.max_episode)));




%%
sim_initial = [0;0;];
noise_power1 = 1;
[local_x_all1, env_x_all1, rect_x_all1, y_xhat_w_v_all1, rl_u_all1] = train.sim([], sim_initial, noise_power1, noise_seed);
[local_x_all2, env_x_all2, rect_x_all2, y_xhat_w_v_all2] = train.sim_lqrcontroller(sim_initial, noise_power1, noise_seed);
    
figure
subplot(2,1,1)
plot(train.t, local_x_all1(:,1),'r-')
hold on,grid on
plot(train.t, local_x_all2(:,1),'r:')
ylabel('\thetta[rad]')
subplot(2,1,2)
plot(train.t, local_x_all1(:,2),'b-')
hold on,grid on
plot(train.t, local_x_all2(:,2),'b:')
ylabel('\omega[frq]')
legend('RL+LQR','LQR')

disp('around equilium phase')
norm(local_x_all1(:,1))
norm(local_x_all2(:,1))
disp('around equilium omega')
norm(local_x_all1(:,2))
norm(local_x_all2(:,2))

%%
noise_seed = 10;

sim_initial = initial;
noise_power1 = 0;
[local_x_all1, env_x_all1, rect_x_all1, y_xhat_w_v_all1, rl_u_all1] = train.sim([], sim_initial, noise_power1, noise_seed);
[local_x_all2, env_x_all2, rect_x_all2, y_xhat_w_v_all2, reward] = train.sim_lqrcontroller(sim_initial, noise_power1, noise_seed);
    
figure
subplot(2,1,1)
plot(train.t, local_x_all1(:,1),'r-')
hold on,grid on
plot(train.t, local_x_all2(:,1),'r:')
ylabel('\theta[rad]')
subplot(2,1,2)
plot(train.t, local_x_all1(:,2),'b-')
hold on,grid on
plot(train.t, local_x_all2(:,2),'b:')
ylabel('\omega[frq]')
legend('RL+LQR','LQR')

disp('train condition phase')
norm(local_x_all1(:,1))
norm(local_x_all2(:,1))
disp('train condition omega')
norm(local_x_all1(:,2))
norm(local_x_all2(:,2))

figure
plot(reward_history)

disp('LQR reward')
reward
