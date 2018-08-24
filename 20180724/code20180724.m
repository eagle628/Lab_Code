% Environment Identification
%{
Nodeの追加したときの雑音に対する影響
retro制御器を新たに増設したとき
シミュレーション中に変化したと仮定
Local＆Envにノイズが両方入っている状況下を考える．
Envノイズの付加される場所にコントローラをを増設する．もしくは外す
今回の実験では，Extendであれば，Ideal Extendを指す
%}
%% initiallize workspace
clear 
close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 30;
seed = 10;
m = [1,2];
d = [2,10]*1e-2;
b = 1;
Y = [1,5];
%Y = [0,1];
r_c = 0.1;
n_ori = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
%% Controlled node number
control_node_o = 2;
%% Identification SN
id_SN = 0.1;
id_noise_node = [2];      % identification noise point
%% Identification Method
id_method = 'ARMAX'; % ARMAX or OE
%% Add node Config
number_add = 1; % number of add node
addnode_initial = [0,0.5];
%% Add Controller COnfig
rng();
addcon_initial = rand(1,6);
%% Add Controller point
control_node_a = 4;
%% Number of Iteration
ITR = 1; 
%% simulation noise point & power
noise_point = control_node_a;
noise_power = 0.1;
%% Save Directory
location = 'C:\Users\Naoya Inoue\Desktop\Test\add_node&controller\noise_E';
%% Node & Figure Name
% Named nodes
ob_v_p  = {strcat('v_node',num2str(control_node_o))};
ob_w_p  = {strcat('w_node',num2str(control_node_o))};
ob_y_p  = {strcat('y_node',num2str(control_node_o))};
ob_y_p2  = {strcat('y_node',num2str(control_node_a))};
local_noise_p = {strcat('d_node',num2str(control_node_o))};
ob_xhat_p = {'xhat_controlled1'}; % Controllerの数がつくべきもの
con_x_p = {'controller_x1'}; % Controllerの数のうちExtendのもののみがつくべきもの
id_noise_p = cell(1,numel(id_noise_node));
for i = 1 : numel(id_noise_node)
    id_noise_p(i) = {strcat('d_node',num2str(id_noise_node(i)))};
end
sim_noise_p = cell(1,numel(noise_point));
for i = 1 :  numel(noise_point)
    sim_noise_p(i) ={strcat('d_node',num2str(noise_point(i)))};
end
% Named experiment condition
%{
name = strcat('N',num2str(Node_number),'_C1node',num2str(control_node_o),'_C2node',num2str(control_node_a),...
                '_addnode_ini'
                );
%}
%% add I/O port for original system
sys_ori = n_ori.get_sys();
add_io = unique([control_node_o,control_node_a,noise_point]);
for i =  1 : numel(add_io)
    sys_ori = n_ori.add_io(sys_ori,add_io(i), strcat('node',num2str(add_io(i))));
end
%% Generate Add node System 
% Generate network for Add node
n_add = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_add.Adj_ref = n_add.Adj_ref*0;
% Make New node
for i = 1 : number_add
    n_add.add_node(m, d, b, Y, r_c);
end
% New Network add I/O port
sys_add = n_add.get_sys();
for i =  1 : numel(add_io)
    sys_add = n_add.add_io(sys_add,add_io(i), strcat('node',num2str(add_io(i))));
end
%% Environment Character (Original)
% Chose Range Fruquency
min = -6;
max =  6;
[~, sys_env_h] = n_ori.get_sys_local(control_node_o);
sys_env_o1 = balred(sys_env_h,6);
[mag_env,pha_env,omega] = bode(sys_env_h,{10^min,10^max});
mag_env=squeeze(mag_env(1,1,:));
pha_env=squeeze(pha_env(1,1,:));
% Environment for add 2 controller
[~, sys_env_o2] = n_ori.get_sys_local(control_node_a);
sys_env_o2 = balred(sys_env_o2,6);
%% Environment Character (add_system)
[~, sys_env_a_h] = n_add.get_sys_local(control_node_o);
sys_env_a = balred(sys_env_a_h,6);
[mag_env_a,pha_env_a,~] = bode(sys_env_a_h,omega);
mag_env_a = squeeze(mag_env_a(1,1,:));
pha_env_a = squeeze(pha_env_a(1,1,:));
%% LQR Controller Parameter
Q = kron(eye(1*numel(control_node_o)),diag([1,1000]));
R = kron(eye(1*numel(control_node_o)),diag([1e-3]));
%% Add Controller for original system (Ideal Extend)
n_ori.controllers = {};
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
sys_ori_extend1 = n_ori.get_sys_controlled(sys_ori);
%% Add Controller for add node system (Ideal Extend)
% Environment model is assumed original system environment.
% Environment system is Changed.
n_add.controllers = {};
n_add.add_controller( control_node_o, sys_env_o1, Q, R);
sys_add_extend1 = n_add.get_sys_controlled(sys_add);
%% Add 2 Controller for original system (Ideal Extend)
n_ori.controllers = {};
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
n_ori.add_controller( control_node_a, sys_env_o2, Q, R);
%n_ori.add_controller( control_node_a, Q, R);
sys_ori_extend2 = n_ori.get_sys_controlled(sys_ori);
%% Add simulation
simulation_time = 100;%simulation_time
Ts_s = 0.01;
t_s = 0:Ts_s:simulation_time-Ts_s;
t_s = t_s';

sim_noise = zeros(length(t_s),4);
sim_seed = 10;
cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',length(t_s),'NumChannels',4,'RandomStream','mt19937ar with seed','Seed',sim_seed);
noise = cn3()*noise_power;
clear cn3
%sim_noise = noise;
sim_noise(:,4) = noise(:,4);
sim_noise(:,2) = noise(:,2);
%sim_noise(round(length(t_s)/2),2) = 1/Ts_s;
%sim_noise(round(length(t_s)/4),3) = 1/Ts_s;
rng();
add_node_initial = rand(1,2);
rng();
add_controller_initial = rand(1);

[y1,xhat1] = sim_change_system(sys_ori_extend1,sys_add_extend1,ob_y_p,{local_noise_p,sim_noise_p},sim_noise,t_s,add_node_initial,25);
[y2,xhat2] = sim_change_system(sys_ori_extend1,sys_ori_extend2,ob_y_p,{local_noise_p,sim_noise_p},sim_noise,t_s,add_controller_initial,25);
y3 = lsim(sys_ori_extend1(ob_y_p,{local_noise_p,sim_noise_p}),sim_noise,t_s);
xhat3 = lsim(sys_ori_extend1('xhat_controlled1',{local_noise_p,sim_noise_p}),sim_noise,t_s);

fig1 = figure_config.set_figure_retro('Add_Response');
state = 2;
% Not Changes
figure_config.plot_retro2( fig1.axes, t_s,y3(:,state), xhat3(:,state), {'g-','linewidth',3.0});
% Add Node
figure_config.plot_retro2( fig1.axes, t_s,y1(:,state), xhat1(:,state), {'b-','linewidth',1.0});
% Add Controller
figure_config.plot_retro2( fig1.axes, t_s,y2(:,state), xhat2(:,state), {'r-','linewidth',1.0});
%% Local Function