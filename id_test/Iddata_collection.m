%% initiallize workspace
clear 
close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 30;
seed = 10;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;

%% Set I/O node 
node = [1,2,3,4];

%% add I/O port for identificaiton
sys_org = n.get_sys();
for i =  1 : numel(node)
    sys_org = n.add_io(sys_org,node(i), strcat('Id_node',num2str(node(i))));
end

%% Set Identificaiton Input
N = 100000;
Ts = 0.25;
%d = randn(N, 10);
t = (0:N-1)'*Ts;
noise_power = 1;
cn = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',10);
d= cn()*noise_power;

%% Iddata collection

lsim_type = 'foh';
%% ノイズが入らないパターン
% d1
v0 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1'}), d(:,1:2), t,lsim_type);
w0 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1'}), d(:,1:2), t,lsim_type);
%% 同じ大きさのノイズが入るパターん
% d1 d2
v1 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1','d_Id_node2'}), d(:,1:4), t,lsim_type);
w1 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1','d_Id_node2'}), d(:,1:4), t,lsim_type);
% d1 d2
v2 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1','d_Id_node2','d_Id_node3'}), d(:,1:6), t,lsim_type);
w2 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1','d_Id_node2','d_Id_node3'}), d(:,1:6), t,lsim_type);
% d1 d2
v3 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1','d_Id_node2','d_Id_node3','d_Id_node4'}), d(:,1:8), t,lsim_type);
w3 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1','d_Id_node2','d_Id_node3','d_Id_node4'}), d(:,1:8), t,lsim_type);

%% save 
save('iddata_ideal','t','Ts','v0','w0')
save('iddata_plus_1_same','t','Ts','v1','w1')
save('iddata_plus_2_same','t','Ts','v2','w2')
save('iddata_plus_3_same','t','Ts','v3','w3')