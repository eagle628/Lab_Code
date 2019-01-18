clear
close all
%% genereate Network
seed = 10;
Node_number = 4;
net1 = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
% net1 = network_swing_simple(Node_number, 1, [2,10]*1e-2, 1, [0,1], 0.1, seed);
% net1 = network_swing_simple(Node_number, [1,100], [2,10]*1e-1, 1, [1,2], 0.9, seed);
net1.Adj_ref = net1.Adj_ref*0;
%% control & noise Node
c_n = 1;
n_n = [];
%% edge adjustment
tmp_idx = find(net1.Adj(c_n, :)~=0);
tmp_adj = 5 + 5*rand(1, length(tmp_idx));
net1.Adj(c_n ,tmp_idx) = tmp_adj;
net1.Adj(tmp_idx, c_n) = tmp_adj;

% net1.remove_edge([1,3]);
net1.plot()
%% Node Name
ob_v_p  = {strcat('v_node',num2str(c_n))};
ob_w_p  = {strcat('w_node',num2str(c_n))};
ob_y_p  = {strcat('y_node',num2str(c_n))};
ID_in_p = {strcat('d_node',num2str(c_n))};
ob_xhat_p = {'xhat_controlled1'}; % The last number should be Controllers Group number.
Noise = cell(1,numel(n_n));
for i = 1 : numel(n_n)
Noise(i) = {strcat('d_node',num2str(n_n(i)))};
end
%% add I/O port for identificaiton
sys_all = net1.get_sys();
add_io_node = unique([c_n,n_n]);
for nnn = add_io_node
    sys_all = net1.add_io(sys_all, nnn, strcat('node',num2str(nnn)));
end

[sys_local, sys_env] = net1.get_sys_local(c_n);
sys_local_vw = sys_local({'w'},{'v'});
%% env character
[mag_ori,~,wout_ori] = bode(sys_env);
wrange = {wout_ori(1),wout_ori(end)};

Loop = loopsens(sys_env,-sys_local_vw);
Loop.Stable
isstable(sys_all(ob_w_p,ID_in_p))
isstable(sys_all(ob_v_p,ID_in_p))
isstable(sys_all)
%% vw generate
Q = diag([1,1]);
R = 1e-3;

net1.add_controller(c_n, Q, R);
sys_s_retro = net1.get_sys_controlled(sys_all);

net1.controllers = {};
model = balred(sys_env, 4);
net1.add_controller(c_n, model, Q, R);
sys_e_retro = net1.get_sys_controlled(sys_all);

[~,rect_G_s] = G_checkz_d(sys_env, sys_local, ss([],[],[],0));
[~,rect_G_e] = G_checkz_d(sys_env, sys_local, model);

figure
bode(rect_G_s,rect_G_e)
figure
hold on
bode(sys_s_retro(sys_s_retro.OutputGroup.y_node1(2),sys_s_retro.InputGroup.d_node1(1)), 'r', rect_G_s,'b')
bode(sys_e_retro(sys_e_retro.OutputGroup.y_node1(2),sys_e_retro.InputGroup.d_node1(1)), 'r:', rect_G_e,'b:')
figure
bode(sys_env,model)
