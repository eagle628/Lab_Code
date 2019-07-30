close all
clear

%% Release GPU memory
gpuDevice(1);

% % % % % %% model based belief state version (Do not work)
% % % % % 
% % % % % seed1 = 6;
% % % % % Node_number = 30;
% % % % % net = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,2], 0.1, seed1);
% % % % % net.Adj_ref = net.Adj_ref*0;
% % % % % 
% % % % % % local system
% % % % % c_n = 1;
% % % % % % edge adjustment
% % % % % tmp_idx = find(net.Adj(c_n, :)~=0);
% % % % % tmp_adj = 2 + 3*rand(1, length(tmp_idx));
% % % % % net.Adj(c_n ,tmp_idx) = tmp_adj;
% % % % % net.Adj(tmp_idx, c_n) = tmp_adj;
% % % % % 
% % % % % % train condition
% % % % % Ts = 0.01;
% % % % % Te = 4;
% % % % % 
% % % % % % model = swing_network_model(net, c_n, Ts);
% % % % % [sys_local, sys_env] = net.get_sys_local(c_n);
% % % % % apx_nx = 6;%net.N*2 -2;
% % % % % % netwrok identificaiton
% % % % % sys_all = net.get_sys();
% % % % % for nnn = 1 : net.N
% % % % %     sys_all = net.add_io(sys_all, nnn, strcat('node',num2str(nnn)));
% % % % % end
% % % % % 
% % % % % net.add_controller(c_n);
% % % % % sys_con = net.get_sys_controlled(sys_all);
% % % % % N = 1000;
% % % % % t = (0:N-1)'*Ts;
% % % % % ob_v_p  = {strcat('v_node',num2str(c_n))};
% % % % % ob_w_p  = {strcat('w_node',num2str(c_n))};
% % % % % ID_in_p = {strcat('d_node',num2str(c_n))};
% % % % % [sys_v_id, ~, Tv, Tiv] = balreal(c2d(sys_con(ob_v_p, ID_in_p), Ts, 'foh'));
% % % % % [sys_w_id, ~, Tw, Tiw] = balreal(c2d(sys_con(ob_w_p, ID_in_p), Ts, 'foh'));
% % % % % rng(12);
% % % % % d_id = randn(N,2);
% % % % % v_id = lsim(sys_v_id, d_id, t);
% % % % % w_id = lsim(sys_w_id, d_id, t);
% % % % % v = v_id + 0;
% % % % % w = w_id + 0;
% % % % % data = struct();
% % % % % data.Ts = Ts;
% % % % % data.u = w;
% % % % % data.y = v;
% % % % % sys_local_vw = sys_local({'w'},{'v'});
% % % % % model_dim = apx_nx;
% % % % % sys_iv = CLRIVC(data,balreal(-sys_local_vw),[model_dim,model_dim+1,model_dim,model_dim],1,5);
% % % % % % identiificaiton system(Env system)
% % % % % sys_iv_G = balreal(ss(sys_iv.G));
% % % % % sys_oe = balreal(d2c(ss(oe(iddata(v, w, Ts), [model_dim+1, model_dim, 0]))));
% % % % % % sys_arx = balreal(d2c(ss(oe(iddata(v, w, Ts), [model_dim, model_dim+1, 0]))));
% % % % % sys_real = balreal(balred(sys_env,apx_nx));
% % % % % % add retro controller
% % % % % net1.controllers = {};
% % % % % 
% % % % % apx_sys = {ss(0), sys_iv_G, sys_real, sys_oe};
% % % % % % apx_sys = {ss(0),sys_iv_G};
% % % % % name = {'simple','iv','real','oe'};
% % % % % %%
% % % % % for idx = 1 : length(apx_sys)
% % % % % try
% % % % %     mail_message(strcat('RL train start', num2str(idx)))
% % % % %     %% Sngle Test
% % % % %     model = swing_network_model(net, c_n, Ts, apx_sys{idx});
% % % % %     % define basis function
% % % % %     basis_N = 3;
% % % % %     range = [-2, 2];
% % % % %     width = (range(2)-range(1))/(basis_N-1);
% % % % %     m = range(1):width:range(2);
% % % % %     nnn = model.local_nx + model.rect_nx;
% % % % % %     nnn = model.rect_nx;
% % % % %     mu = m;
% % % % %     for itr = 1 : nnn-1
% % % % %         mu = combvec(mu, m); % col vector combinater function
% % % % %     end
% % % % %     mu = mu';
% % % % %     sigma = 1*ones(size(mu, 1), 1);
% % % % %     RBF1 = Radial_Basis_Function(size(mu, 1), mu, sigma);
% % % % %     % RBF1 = Normalized_Radial_Basis_Function(size(mu, 1), mu, sigma);
% % % % %     % policy
% % % % %     sigma_pi = 0.5;
% % % % %     policy = policy_RBF(RBF1, sigma_pi);
% % % % %     % value
% % % % %     value  =  value_RBF(RBF1);
% % % % % 
% % % % %     train_seed = 1024;
% % % % %     train = netwrok_retro_by_actor_critic_with_eligibility_traces_episodic(model, policy, value, Te);
% % % % % 
% % % % %     initial = [0,0.1];
% % % % % %     initial = [0;1;];
% % % % %     [local_x_all, mpc_u_all, rl_u_all, theta_mu_snpashot, theta_sigma_snpashot, w_snpashot, reward_history] = train.train(initial, train_seed);
% % % % %     %%
% % % % % % %     sim_initial = [0;0;];
% % % % % % %     [local_x_all1, env_x_all1, rect_x_all1, y_xhat_w_v_all1, rl_u_all1] = train.sim([], sim_initial, 6);
% % % % % % %     [local_x_all2, env_x_all2, rect_x_all2, y_xhat_w_v_all2] = train.sim_lqrcontroller(sim_initial, 6);
% % % % %     
% % % % % %     save(strcat('RL_net',num2str(Node_number),'_discount99_', name{idx}, '_reward_fix_Q_episode_1000_No2'));
% % % % % %     save(strcat('RL_net',num2str(Node_number),'_discount99_', name{idx}, '_reward_fix_Q_episode_2000_NonLQR'));
% % % % % %     save(strcat('RL_net',num2str(Node_number),'_discount99_', name{idx}, '_reward_only_local_episode_',num2str(train.max_episode)));
% % % % %     save('test')
% % % % % catch ME
% % % % %     mail_message(jsonencode(ME));
% % % % % end
% % % % % mail_message(strcat('RL train End', num2str(idx)));
% % % % % end
% % % % % mail_message('RL train End All');


clear
close all
seed1 = 6;
Node_number = 30;
net = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,2], 0.1, seed1);
net.Adj_ref = net.Adj_ref*0;
% local system
c_n = 1;
% edge adjustment
tmp_idx = find(net.Adj(c_n, :)~=0);
tmp_adj = 2 + 3*rand(1, length(tmp_idx));
net.Adj(c_n ,tmp_idx) = tmp_adj;
net.Adj(tmp_idx, c_n) = tmp_adj;
% net.plot()
% train condition
Ts = 0.01;
Te = 4;
belief_N = 8;
basis_N = 3;
range = [-2,2];
width = (range(2)-range(1))/(basis_N-1);
m = range(1):width:range(2);
nnn = belief_N;
mu = m;
for itr = 1 : nnn-1
mu = combvec(mu, m); % col vector combinater function
end
mu = mu';
sigma = 1*ones(size(mu, 1), 1);
RBF1 = Radial_Basis_Function(size(mu, 1), mu, sigma);
% RBF1 = Normalized_Radial_Basis_Function(size(mu, 1), mu, sigma);
% policy
sigma_pi = 0.5;
policy = policy_RBF(RBF1, sigma_pi);
% value
value  =  value_RBF(RBF1);
model = swing_network_model(net, c_n, Ts);
train_seed = 1024;
train = netwrok_retro_by_actor_critic_with_eligibility_traces_episodic(model, policy, value, Te, belief_N);
initial = [0,0];
[local_x_all, mpc_u_all, rl_u_all, theta_mu_snpashot, theta_sigma_snpashot, w_snpashot, reward_history] = train.train(initial, train_seed, 'parallel', 'off');


%%
noise_seed = 10;

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

sim_initial = [0;1];
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
