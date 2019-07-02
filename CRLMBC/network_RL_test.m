close all
clear

seed1 = 6;
net = network_swing_simple(10, [1,2], [2,10]*1e-2, 1, [1,2], 0.1, seed1);
c_n = 1;
Ts = 0.01;

% model = swing_network_model(net, c_n, Ts);
[sys_local, sys_env] = net.get_sys_local(c_n);
model = swing_network_model(net, c_n, Ts, sys_env);

seed2 = 1024;
train = netwrok_retro_by_actor_critic_with_eligibility_traces_episodic(model, 1, 2, seed2);

train.train([1;0]);
