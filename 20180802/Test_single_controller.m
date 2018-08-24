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
close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 10;
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
% n_ori.plot()
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
control_node_a1 = 2;
control_node_a2 = 3;
%% Number of Iteration
ITR = 1; 
%% simulation noise point & power
noise_point = [control_node_a1,control_node_a2,4];
% noise_power = 0.1;
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
%% Environment Character (Original)
% Chose Range Fruquency
min = -6;
max =  6;
[~, sys_env_o1] = n_ori.get_sys_local(control_node_o);
[mag_env,pha_env,omega] = bode(sys_env_o1,{10^min,10^max});
mag_env=squeeze(mag_env(1,1,:));
pha_env=squeeze(pha_env(1,1,:));
% sys_env_o1 = balred(sys_env_h,6);
% Environment for add 2 controller
[~, sys_env_o2] = n_ori.get_sys_local(control_node_a1);
% sys_env_o2 = balred(sys_env_o2,6);
% Environment for add 3 controller
[~, sys_env_o3] = n_ori.get_sys_local(control_node_a2);
% sys_env_o3 = balred(sys_env_o3,6);
%% LQR Controller Parameter
Q = kron(eye(1*numel(control_node_o)),diag([1,1000]));
R = kron(eye(1*numel(control_node_o)),diag([1e-3]));
%% SS Controller for original system (Ideal Extend)
n_ori.controllers = {};
n_ori.add_controller( control_node_o, Q, R);
n_ori.add_controller( control_node_a1, Q, R);
sys_ori_c2_ss = n_ori.get_sys_controlled(sys_ori);
%% EE Controller for original system (Ideal Extend)
n_ori.controllers = {};
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
n_ori.add_controller( control_node_a1,sys_env_o1, Q, R);
sys_ori_c2_ee = n_ori.get_sys_controlled(sys_ori);
%% ES Controller for original system (Ideal Extend)
n_ori.controllers = {};
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
n_ori.add_controller( control_node_a1, Q, R);
sys_ori_c2_es = n_ori.get_sys_controlled(sys_ori);

%% Add simulation
simulation_time = 100;%simulation_time 
ch_t = 100;
Ts_s = 0.01;
t_s = 0:Ts_s:simulation_time-Ts_s;
t_s = t_s';

sim_noise = zeros(length(t_s),8);
rng('shuffle');
%sim_seed = 10;
sim_seed = randi(1000,1,1);

cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',length(t_s),'NumChannels',8,'RandomStream','mt19937ar with seed','Seed',sim_seed);
noise = cn3();
clear cn3
%sim_noise = noise;
noise_power = 0;
sim_noise(:,8) = noise(:,8)*noise_power;
sim_noise(:,6) = noise(:,6)*noise_power;
sim_noise(:,4) = noise(:,4)*noise_power;
sim_noise(:,2) = noise(:,2);

Noise_port = {local_noise_p{:},sim_noise_p{:}};

% SS
y_ss = lsim(sys_ori_c2_ss(ob_y_p,Noise_port),sim_noise,t_s,'zoh');
xhat_ss = lsim(sys_ori_c2_ss(ob_xhat_p,Noise_port),sim_noise,t_s);
v_ss = lsim(sys_ori_c2_ss(ob_v_p,Noise_port),sim_noise,t_s);
vhat_con_ss = lsim(sys_ori_c2_ss(ob_vhat_p,Noise_port),sim_noise,t_s);
w_ss = lsim(sys_ori_c2_ss(ob_w_p,Noise_port),sim_noise,t_s);
vhat_ss = lsim( sys_env_o1, w_ss, t_s);
% EE
y_ee = lsim(sys_ori_c2_ee(ob_y_p,Noise_port),sim_noise,t_s,'zoh');
xhat_ee = lsim(sys_ori_c2_ee(ob_xhat_p,Noise_port),sim_noise,t_s);
v_ee = lsim(sys_ori_c2_ee(ob_v_p,Noise_port),sim_noise,t_s);
vhat_con_ee = lsim(sys_ori_c2_ee(ob_vhat_p,Noise_port),sim_noise,t_s);
w_ee = lsim(sys_ori_c2_ee(ob_w_p,Noise_port),sim_noise,t_s);
vhat_ee = lsim( sys_env_o1, w_ee, t_s);
% ES
y_es = lsim(sys_ori_c2_ee(ob_y_p,Noise_port),sim_noise,t_s,'zoh');
xhat_es = lsim(sys_ori_c2_ee(ob_xhat_p,Noise_port),sim_noise,t_s);
v_es = lsim(sys_ori_c2_ee(ob_v_p,Noise_port),sim_noise,t_s);
vhat_con_es = lsim(sys_ori_c2_ee(ob_vhat_p,Noise_port),sim_noise,t_s);
w_es = lsim(sys_ori_c2_ee(ob_w_p,Noise_port),sim_noise,t_s);
vhat_es = lsim( sys_env_o1, w_es, t_s);
%% InterConnection Signal
fig2 = figure_config.set_figure_retro('InterConnection v',[0,0,900,900],{'v','\hat{v}','\tilde{v}'});
% Not Changes
vtilde_ss = figure_config.plot_retro1( fig2.axes, t_s, v_ss, vhat_ss, {'g-','linewidth',3.0});
% Add Node
vtilde_ee = figure_config.plot_retro1( fig2.axes, t_s, v_ee, vhat_ee, {'b-','linewidth',1.5});
% Add Controller
vtilde_es = figure_config.plot_retro1( fig2.axes, t_s, v_es, vhat_es, {'r-','linewidth',0.8});

SS = fig2.axes.ax1.Children(end);
EE = fig2.axes.ax1.Children(end-1);
ES = fig2.axes.ax1.Children(end-2);
%ex2_int = fig2.axes.ax1.Children(end-4);
legend(fig2.axes.ax1,[SS,EE,ES],...
    'SS','EE','ES',...
    'location','best')
