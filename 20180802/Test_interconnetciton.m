% Environment Identification
%{
Nodeの追加したときの雑音に対する影響
retro制御器を新たに増設したとき
シミュレーション中に変化したと仮定
Local＆Envにノイズが両方入っている状況下を考える．
Envノイズの付加される場所にコントローラをを増設する．もしくは外す
今回の実験では，Extendであれば，Ideal Extendを指す
%%%% 今更だが，xhatは，結局無視する量を示すので，ytildeと同値になるはず，正確には，Cがeye(1)だからだが，

v \tilde{v} \hat{v}との相関を考える（時間区切りPEMを目指す）
%}
%% initiallize workspace
clear 
% close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 4;
seed = 10;
%rng('shuffle');
%seed = randi(1000,1,1);
m = [1,2];
d = [2,10]*1e-2;
b = 1;
Y = [1,5];
% Y = [0,1];
r_c = 0.1;
n_ori = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% Controlled node number
control_node_o = 1;
%% Identification SN
id_SN = 0.1;
id_noise_node = [2];      % identification noise point
%% Identification Method
id_method = 'ARMAX'; % ARMAX or OE
%% Add node Number
number_add = 1; % number of add node
%% Add Controller point
control_node_a1 = 4;
control_node_a2 = 3;
%% Number of Iteration
ITR = 1; 
%% simulation noise point & power
noise_point = [control_node_a1,control_node_a2,1];
noise_power = 0.1;
%% Save Directory
location = 'C:\Users\Naoya Inoue\Desktop\Test\add_node&controller\noise_E';
%% Node & Figure Name
% Named nodes
ob_v_p  = {strcat('v_node',num2str(control_node_o))};
ob_w_p  = {strcat('w_node',num2str(control_node_o))};
ob_y_p  = {strcat('y_node',num2str(control_node_o))};
ob_y_p2  = {strcat('y_node',num2str(control_node_a1))};
local_noise_p = {strcat('d_node',num2str(control_node_o))};
ob_xhat_p = {'xhat_controlled1'}; % Controllerの数がつくべきもの
ob_vhat_p = {'vhat_controlled1'}; % Controllerの数がつくべきもの
con_x_p = {'controller_x1'}; % Controllerの数のうちExtendのもののみがつくべきもの
% identificaiton noise point
id_noise_p = cell(1,numel(id_noise_node));
for i = 1 : numel(id_noise_node)
    id_noise_p(i) = {strcat('d_node',num2str(id_noise_node(i)))};
end
% simlation noise point
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
add_io = unique([control_node_o,control_node_a1,control_node_a2,noise_point]);
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
n_add.plot()
%% add I/O port for Add node system
sys_add = n_add.get_sys();
for i =  1 : numel(add_io)
    sys_add = n_add.add_io(sys_add,add_io(i), strcat('node',num2str(add_io(i))));
end
%% Environment Character (Original)
% Chose Range Fruquency
min = -6;
max =  6;
[~, sys_env_h] = n_ori.get_sys_local(control_node_o);
% sys_env_o1 = balred(sys_env_h,6);
sys_env_o1 = sys_env_h;
[mag_env,pha_env,omega] = bode(sys_env_h,{10^min,10^max});
mag_env=squeeze(mag_env(1,1,:));
pha_env=squeeze(pha_env(1,1,:));
%
[~, sys_env_h] = n_ori.get_sys_local(control_node_o);
% sys_env_o1 = balred(sys_env_h,6);
sys_env_o2 = sys_env_h;
%% LQR Controller Parameter
Q = kron(eye(1*numel(control_node_o)),diag([1,1000]));
R = kron(eye(1*numel(control_node_o)),diag([1e-3]));
%% Add 1 Controller for original system (Ideal Extend)
n_ori.controllers = {};
% n_ori.add_controller( control_node_o, Q, R);
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
sys_ori_extend1 = n_ori.get_sys_controlled(sys_ori);
%% add 2 controller for original system
n_ori.controllers = {};
% n_ori.add_controller( control_node_o, Q, R);
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
% n_ori.add_controller( 2, Q, R);
n_ori.add_controller( 2, sys_env_o2, Q, R);
sys_ori_extend2 = n_ori.get_sys_controlled(sys_ori);
%% add 1 controller for add node system
n_add.controllers = {};
% n_add.add_controller( control_node_o, Q, R);
n_add.add_controller( control_node_o, sys_env_o1, Q, R);
sys_add_extend1 = n_add.get_sys_controlled(sys_add);
%% simulation
simulation_time =1000;%simulation_time 
ch_t = 100;
Ts_s = 0.01;
t_s = 0:Ts_s:simulation_time-Ts_s;
t_s = t_s';

Noise_port = {local_noise_p{:},sim_noise_p{:}};
sim_noise = zeros(length(t_s),numel(Noise_port)*2);
rng('shuffle');
sim_seed = 10;
% sim_seed = randi(1000,1,1);

cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',length(t_s),'NumChannels',numel(noise_point)*8,'RandomStream','mt19937ar with seed','Seed',sim_seed);
noise = cn3()*noise_power;
clear cn3

% sim_noise(:,8) = noise(:,8);
% sim_noise(:,6) = noise(:,6);
% sim_noise(:,4) = noise(:,4);
sim_noise(:,2) = noise(:,2);
%%%%%%%%%%%%%%%%%
rng('shuffle');
% add_controller_initial = (rand(1,6)-0.5);
add_controller_initial = ones(1,6);
rng('shuffle');
% add_node_initial = (rand(1,2)-0.5);
add_node_initial = ones(1,2);
% add controller
[y1,xhat1,v1,vhat_con1,w1] = sim_change_system(1,sys_ori_extend1,sys_ori_extend2,ob_y_p,Noise_port,sim_noise,t_s,add_controller_initial,ch_t);
vhat1 = lsim( sys_env_o1, w1, t_s);
% add node
[y2,xhat2,v2,vhat_con2,w2] = sim_change_system(1,sys_ori_extend1,sys_add_extend1,ob_y_p,Noise_port,sim_noise,t_s,add_node_initial,ch_t);
vhat2 = lsim( sys_env_o1, w2, t_s);


fig1 = figure_config.set_figure_retro('InterConnection v',[0,0,900,900],{'v','\hat{v}','\tilde{v}'});
% Not Changes
vtilde1 = figure_config.plot_retro1( fig1.axes, t_s, v1, vhat1, {'r-','linewidth',0.8});
vtilde2 = figure_config.plot_retro1( fig1.axes, t_s, v2, vhat2, {'b-','linewidth',0.8});

legend(fig1.axes.ax1, [fig1.axes.ax1.Children(2),fig1.axes.ax1.Children(1)],'add controller','add node')
% figure
% subplot(3,1,1)
% plot(t_s,v1)
% subplot(3,1,2)
% plot(t_s,vhat1)
% subplot(3,1,3)
% plot(t_s,v1-vhat1)