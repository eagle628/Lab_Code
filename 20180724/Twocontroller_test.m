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
Node_number = 4;
seed = 10;
m = [1,2];
d = [2,10]*1e-2;
b = 1;
Y = [1,5];
%Y = [0,1];
r_c = 0.1;
n_ori = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
%% Add node Config
number_add = 1; % number of add node
add_initial = [1,0];
%% Generate Add node System 
% Generate network for Add node
n_add = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_add.Adj_ref = n_add.Adj_ref*0;
% Make New node
for i = 1 : number_add
    n_add.add_node(m, d, b, Y, r_c);
end
%% Controlled node number
control_node_o = 1;
%% Add Controller point
control_node_a = 4;
%% Identification SN
id_SN = 0.1;
id_noise_node = 3;
%% Identification Method
id_method = 'ARMAX'; % ARMAX or OE
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
name_c = '';
for i = 1 : numel(control_node_o)
    name_c = strcat(name_c,num2str(control_node_o(i)),'+');
end
name_n = '';
for i = 1 : numel(id_noise_node)
    name_n = strcat(name_n,num2str(id_noise_node(i)),'+');
end
str1 = replace(num2str(id_SN),'.','%');
str2 = replace(num2str(noise_power),'.','-');
%{
name = strcat('N',num2str(Node_number),'_Cnode',num2str(control_node_o),'_ACnode',num2str(control_node_a),...
                '_NoisePower',str2,'_NoiseLocation_',num2str(noise_point)...
                );
%}
%% add I/O port for original system
sys_ori = n_ori.get_sys();
add_io = unique([control_node_o,control_node_a,noise_point]);
for i =  1 : numel(add_io)
    sys_ori = n_ori.add_io(sys_ori,add_io(i), strcat('node',num2str(add_io(i))));
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
n_ori.controllers = [];
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
%n_ori.add_controller( control_node_o, Q, R);
sys_ori_extend1 = n_ori.get_sys_controlled(sys_ori);
%% Add Controller for original system (Ideal Extend)
n_ori.controllers = [];
n_ori.add_controller( control_node_a, sys_env_o2, Q, R);
%n_ori.add_controller( control_node_o, Q, R);
sys_ori_extend1_a = n_ori.get_sys_controlled(sys_ori);
%% Add 2 Controller for original system (Ideal Extend)
n_ori.controllers = [];
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
n_ori.add_controller( control_node_a, sys_env_o2, Q, R);
%n_ori.add_controller( control_node_o, Q, R);
%n_ori.add_controller( control_node_a, Q, R);
sys_ori_extend2 = n_ori.get_sys_controlled(sys_ori);
%%
simulation_time = 100;%simulation_time
Ts_s = 0.01;
t_s = 0:Ts_s:simulation_time-Ts_s;
t_s = t_s';
%%%%%%%%
figure('Name','2 controller','position',[-1800,700,600,300])
y = impulse(sys_ori_extend2(ob_y_p,local_noise_p),t_s);
plot(t_s,y(:,2,1),':','linewidth',2.0)
hold on
y = impulse(sys_ori_extend2(ob_y_p2,sim_noise_p),t_s);
plot(t_s,y(:,2,1))
%%%%%%%%
figure('Name','1 controller_o','position',[-1200,700,600,300])
y = impulse(sys_ori_extend1(ob_y_p,local_noise_p),t_s);
plot(t_s,y(:,2,1),':','linewidth',2.0)
hold on
y = impulse(sys_ori_extend1(ob_y_p2,sim_noise_p),t_s);
plot(t_s,y(:,2,1))
%%%%%%%%
figure('Name','1 controller_a','position',[-1200,300,600,300])
y = impulse(sys_ori_extend1_a(ob_y_p,local_noise_p),t_s);
plot(t_s,y(:,2,1),':','linewidth',2.0)
hold on
y = impulse(sys_ori_extend1_a(ob_y_p2,sim_noise_p),t_s);
plot(t_s,y(:,2,1))
%%%%%%%%
figure('Name','Original','position',[-600,700,600,300])
y = impulse(sys_ori(ob_y_p,local_noise_p),t_s);
plot(t_s,y(:,2,1),':','linewidth',2.0)
hold on
y = impulse(sys_ori(ob_y_p2,sim_noise_p),t_s);
plot(t_s,y(:,2,1))