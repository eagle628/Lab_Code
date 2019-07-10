close all
clear

%% Release GPU memory
gpuDevice(1)

%%

seed1 = 6;
Node_number = 4;
net = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,2], 0.1, seed1);
c_n = 1;
Ts = 0.01;
Te = 10;

% model = swing_network_model(net, c_n, Ts);
[sys_local, sys_env] = net.get_sys_local(c_n);
apx_nx = 6;
% netwrok identificaiton
sys_all = net.get_sys();
for nnn = 1 : net.N
    sys_all = net.add_io(sys_all, nnn, strcat('node',num2str(nnn)));
end

net.add_controller(c_n);
sys_con = net.get_sys_controlled(sys_all);
N = 1000;
t = (0:N-1)'*Ts;
ob_v_p  = {strcat('v_node',num2str(c_n))};
ob_w_p  = {strcat('w_node',num2str(c_n))};
ID_in_p = {strcat('d_node',num2str(c_n))};
[sys_v_id, ~, Tv, Tiv] = balreal(c2d(sys_con(ob_v_p, ID_in_p), Ts, 'foh'));
[sys_w_id, ~, Tw, Tiw] = balreal(c2d(sys_con(ob_w_p, ID_in_p), Ts, 'foh'));
d_id = randn(N,2);
v_id = lsim(sys_v_id, d_id, t);
w_id = lsim(sys_w_id, d_id, t);
v = v_id + 0;
w = w_id + 0;
data = struct();
data.Ts = Ts;
data.u = w;
data.y = v;
sys_local_vw = sys_local({'w'},{'v'});
model_dim = apx_nx;
sys = CLRIVC(data,balreal(-sys_local_vw),[model_dim,model_dim+1,model_dim,model_dim],1,5);

net1.add_controller = {};

%% RL

model = swing_network_model(net, c_n, Ts, balreal(balred(sys_env,apx_nx)));

seed2 = 1024;
basis_N = 3;
train = netwrok_retro_by_actor_critic_with_eligibility_traces_episodic(model, Te, basis_N, seed2);

[local_x_all, mpc_u_all, rl_u_all, theta_mu, w, reward_history] = train.train([1;0]);
%%
[local_x_all1, env_x_all1, rect_x_all1, y_xhat_w_v_all1, rl_u_all1] = train.sim([], 6);
[local_x_all2, env_x_all2, rect_x_all2, y_xhat_w_v_all2] = train.sim_lqrcontroller(6);

figure
subplot(2,1,1)
plot(train.t, local_x_all1(:,1),'r-')
hold on
plot(train.t, local_x_all2(:,1),'r:')
subplot(2,1,2)
plot(train.t, local_x_all1(:,2),'b-')
hold on
plot(train.t, local_x_all2(:,2),'b:')
