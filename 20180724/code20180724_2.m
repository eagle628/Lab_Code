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
%function [y1,y2,y3,y4,xhat1,xhat2,xhat3,xhat4,t_s] = code20180724_2()
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 4;
seed = 8;
%rng('shuffle');
%seed = randi(1000,1,1);
m = [1,2];
d = [2,10]*1e-2;
b = 1;
Y = [1,5];
%Y = [0,1];
r_c = 0.1;
n_ori = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% Controlled node number
control_node_o = 2;
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
noise_point = [control_node_a1,control_node_a2];
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
sys_env_o1 = balred(sys_env_h,6);
[mag_env,pha_env,omega] = bode(sys_env_h,{10^min,10^max});
mag_env=squeeze(mag_env(1,1,:));
pha_env=squeeze(pha_env(1,1,:));
% Environment for add 2 controller
[~, sys_env_o2] = n_ori.get_sys_local(control_node_a1);
sys_env_o2 = balred(sys_env_o2,6);
% Environment for add 3 controller
[~, sys_env_o3] = n_ori.get_sys_local(control_node_a2);
sys_env_o3 = balred(sys_env_o3,6);
%% Environment Character (add_system)
[~, sys_env_a_h] = n_add.get_sys_local(control_node_o);
sys_env_a = balred(sys_env_a_h,6);
[mag_env_a,pha_env_a,~] = bode(sys_env_a_h,omega);
mag_env_a = squeeze(mag_env_a(1,1,:));
pha_env_a = squeeze(pha_env_a(1,1,:));
%% LQR Controller Parameter
Q = kron(eye(1*numel(control_node_o)),diag([1,1000]));
R = kron(eye(1*numel(control_node_o)),diag([1e-3]));
%% Add 2 Controller for original system (Ideal Extend)
n_ori.controllers = {};
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
n_ori.add_controller( control_node_a1, sys_env_o2, Q, R);
sys_ori_extend2 = n_ori.get_sys_controlled(sys_ori);
%% Add 2 Controller for add system (Ideal Extend)
n_add.controllers = {};
n_add.add_controller( control_node_o, sys_env_o1, Q, R);
n_add.add_controller( control_node_a1, sys_env_o2, Q, R);
%n_add.add_controller( Node_number+1,  Q, R);
sys_add_extend2 = n_add.get_sys_controlled(sys_add);
%% Add 3 Controller for original system (Ideal Extend)
n_ori.controllers = {};
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
n_ori.add_controller( control_node_a1, sys_env_o2, Q, R);
n_ori.add_controller( control_node_a2, sys_env_o3, Q, R);
sys_ori_extend3 = n_ori.get_sys_controlled(sys_ori);
%% Remove Controller for original system
% 2番目のコントローラを排除
n_ori.controllers = {};
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);
sys_ori_extend1 = n_ori.get_sys_controlled(sys_ori);
%% Add simulation
simulation_time =1000;%simulation_time 
ch_t = 25;
Ts_s = 0.01;
t_s = 0:Ts_s:simulation_time-Ts_s;
t_s = t_s';

sim_noise = zeros(length(t_s),4);
rng('shuffle');
%sim_seed = 10;
sim_seed = randi(1000,1,1);

cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',length(t_s),'NumChannels',6,'RandomStream','mt19937ar with seed','Seed',sim_seed);
noise = cn3()*noise_power;
clear cn3
%sim_noise = noise;
sim_noise(:,6) = noise(:,6);
sim_noise(:,4) = noise(:,4);
sim_noise(:,2) = noise(:,2);
%sim_noise(round(length(t_s)/2),2) = 1/Ts_s;
%sim_noise(round(length(t_s)/4),3) = 1/Ts_s;
rng('shuffle');
add_node_initial = (rand(1,2)-0.5)/5;
%add_node_initial = [0,0];
rng('shuffle');
add_controller_initial = (rand(1,6)-0.5)/5;
%add_controller_initial = 0;

Noise_port = {local_noise_p{:},sim_noise_p{:}};

[y1,xhat1] = sim_change_system(sys_ori_extend2,sys_add_extend2,ob_y_p,Noise_port,sim_noise,t_s,add_node_initial,ch_t);
[y2,xhat2] = sim_change_system(sys_ori_extend2,sys_ori_extend3,ob_y_p,Noise_port,sim_noise,t_s,add_controller_initial,ch_t);
y3 = lsim(sys_ori_extend2(ob_y_p,Noise_port),sim_noise,t_s);
xhat3 = lsim(sys_ori_extend2('xhat_controlled1',Noise_port),sim_noise,t_s);
[y4,xhat4] = sim_change_system(sys_ori_extend2,sys_ori_extend1,ob_y_p,Noise_port,sim_noise,t_s,[2],25);
ch_t_p = round(ch_t/Ts_s)+1;
sim_noise(ch_t_p,2) = 1/Ts_s;
y5 = lsim(sys_ori_extend2(ob_y_p,Noise_port),sim_noise,t_s);
xhat5 = lsim(sys_ori_extend2('xhat_controlled1',Noise_port),sim_noise,t_s);
%
fig1 = figure_config.set_figure_retro('Response',[0,0,900,900]);
state = 2;
% Not Changes
figure_config.plot_retro2( fig1.axes, t_s,y3(:,state), xhat3(:,state), {'g-','linewidth',3.0});
% Add Node
figure_config.plot_retro2( fig1.axes, t_s,y1(:,state), xhat1(:,state), {'b-','linewidth',0.8});
% Add Controller
figure_config.plot_retro2( fig1.axes, t_s,y2(:,state), xhat2(:,state), {'r-','linewidth',0.8});
% Remove Controller
figure_config.plot_retro2( fig1.axes, t_s,y4(:,state), xhat4(:,state), {'k-','linewidth',0.8});
% Internal Noise
%figure_config.plot_retro2( fig1.axes, t_s,y5(:,state), xhat5(:,state), {'c-','linewidth',0.8});

ex2_ori = fig1.axes.ax1.Children(end);
ex2_add = fig1.axes.ax1.Children(end-1);
ex3_ori = fig1.axes.ax1.Children(end-2);
ex1_ori = fig1.axes.ax1.Children(end-3);
%ex2_int = fig1.axes.ax1.Children(end-4);
legend(fig1.axes.ax1,[ex2_ori,ex2_add,ex3_ori,ex1_ori],...
    'Original(Extend2)','AddNode(Extend2)','Original(Extend3)','Original(Extend1)',...
    'location','best')
%}
%end
%% FFT
Fs = 1/Ts_s;
ch_t = ch_t/Ts_s+1;
[ P1_1, P2_1, L] = fft_func( y1(ch_t:length(y1),state), 3);
[ P1_2, P2_2] = fft_func( y2(ch_t:length(y2),state), 3);
[ P1_3, P2_3] = fft_func( y3(ch_t:length(y3),state), 3);
[ P1_4, P2_4] = fft_func( y4(ch_t:length(y4),state), 3);

%[ P1_1, P2_1, L] = fft_func( xhat1(ch_t:length(y1),state), 3);
%[ P1_2, P2_2] = fft_func( xhat2(ch_t:length(y2),state), 3);
%[ P1_3, P2_3] = fft_func( xhat3(ch_t:length(y3),state), 3);
%[ P1_4, P2_4] = fft_func( xhat4(ch_t:length(y4),state), 3);

figure('Name','Power Spectram')
f = Fs*(0:(L/2))/L;
plot(f,Normalization(P1_3),'g-','linewidth',3.0)
hold on
grid on
plot(f,normalize(P1_1,'range'),'b-','linewidth',1.0)
plot(f,normalize(P1_2,'range'),'r-','linewidth',1.0)
plot(f,normalize(P1_4,'range'),'k-','linewidth',1.0)
title('Single-Sided Amplitude Spectrum')
xlabel('f (Hz)')
xlim([0 1])
legend('Original(Extend2)','AddNode(Extend2)','Original(Extend3)','Original(Extend1)')