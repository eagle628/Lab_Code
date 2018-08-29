clear 
close all
%% Generate Network
seed = 3;
Node_number = 2;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;
% n.plot()
%% control & noise Node
c_n = 1;
n_n = [2];
%% signal power
id_in_p = 1;
noise_p = 1;
%% simlationi noise seed
sim_seed = 10;
%% init system method
id_method = 'ARX';
%% Node Name
ob_v  = {strcat('v_node',num2str(c_n))};
ob_w  = {strcat('w_node',num2str(c_n))};
ID_in = {strcat('d_node',num2str(c_n))};
Noise = cell(1,numel(n_n));
for i = 1 : numel(n_n)
Noise(i) = {strcat('d_node',num2str(n_n(i)))};
end
%% add I/O port for identificaiton
sys_org = n.get_sys();
for i =  1 : n.N
sys_org = n.add_io(sys_org,i, strcat('node',num2str(i)));
end
[sys_local, sys_env] = n.get_sys_local(c_n);

all = loopsens(sys_env,-sys_local({'w'},{'v'}));

eig(sys_org.A)
all.Poles

G = sys_env*sys_local({'w'},{'v'})/(1-sys_env*sys_local({'w'},{'v'}));
G = minreal(G);
eig(G.A)