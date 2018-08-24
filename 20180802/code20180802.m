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
%% Generate Add node System 
% Generate network for Add node
n_add = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_add.Adj_ref = n_add.Adj_ref*0;
% Make New node
for i = 1 : number_add
    n_add.add_node(m, d, b, Y, r_c);
end
% n_add.plot()
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
% n_ori.add_controller( control_node_a1, Q, R);
sys_ori_extend2 = n_ori.get_sys_controlled(sys_ori);
%% Add 2 Controller for add system (Ideal Extend)
n_add.controllers = {};
n_add.add_controller( control_node_o, sys_env_o1, Q, R);

n_add.add_controller( control_node_a1, sys_env_o2, Q, R);
% n_add.add_controller( control_node_a1,  Q, R);

sys_add_extend2 = n_add.get_sys_controlled(sys_add);
%% Add 3 Controller for original system (Ideal Extend)
n_ori.controllers = {};
n_ori.add_controller( control_node_o, sys_env_o1, Q, R);

n_ori.add_controller( control_node_a1, sys_env_o2, Q, R);
% n_ori.add_controller( control_node_a1, Q, R);
n_ori.add_controller( control_node_a2, sys_env_o3, Q, R);
% n_ori.add_controller( control_node_a2, Q, R);
sys_ori_extend3 = n_ori.get_sys_controlled(sys_ori);
%% Remove Controller for original system
% 2番目のコントローラを排除
n_ori.controllers = {};
 n_ori.add_controller( control_node_o, sys_env_o1, Q, R);

sys_ori_extend1 = n_ori.get_sys_controlled(sys_ori);
%% Add simulation
simulation_time = 1000;%simulation_time 
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
noise_power = 1;
sim_noise(:,8) = noise(:,8)*noise_power;
sim_noise(:,6) = noise(:,6)*noise_power;
sim_noise(:,4) = noise(:,4)*noise_power;
sim_noise(:,2) = noise(:,2);
% sim_noise = noise;
% sim_noise(:,3:end) = sim_noise(:,3:end)*noise_power;

rng('shuffle');
add_node_initial = (rand(1,2)-0.5);
%add_node_initial = [0,0];
rng('shuffle');
add_controller_initial = (rand(1,6)-0.5);
%add_controller_initial = 0;

Noise_port = {local_noise_p{:},sim_noise_p{:}};

% add node system
[y1,xhat1,v1,vhat_con1,w1] = sim_change_system(1,sys_ori_extend2,sys_add_extend2,ob_y_p,Noise_port,sim_noise,t_s,add_node_initial,ch_t);
% add controller
[y2,xhat2,v2,vhat_con2,w2] = sim_change_system(1,sys_ori_extend2,sys_ori_extend3,ob_y_p,Noise_port,sim_noise,t_s,add_controller_initial,ch_t);
% Not Change
y3 = lsim(sys_ori_extend2(ob_y_p,Noise_port),sim_noise,t_s);
xhat3 = lsim(sys_ori_extend2(ob_xhat_p,Noise_port),sim_noise,t_s);
v3 = lsim(sys_ori_extend2(ob_v_p,Noise_port),sim_noise,t_s);
vhat_con3 = lsim(sys_ori_extend2(ob_vhat_p,Noise_port),sim_noise,t_s);
w3 = lsim(sys_ori_extend2(ob_w_p,Noise_port),sim_noise,t_s);
% Remove Controller
[y4,xhat4,v4,vhat_con4,w4] = sim_change_system(1,sys_ori_extend2,sys_ori_extend1,ob_y_p,Noise_port,sim_noise,t_s,[2],ch_t);

ch_t_p = round(ch_t/Ts_s)+1;
sim_noise(ch_t_p,2) = 1/Ts_s;
% internal Noise
y5 = lsim(sys_ori_extend2(ob_y_p,Noise_port),sim_noise,t_s);
xhat5 = lsim(sys_ori_extend2('xhat_controlled1',Noise_port),sim_noise,t_s);
%
fig1 = figure_config.set_figure_retro('Response',[0,0,900,900]);
state = 2;
% Not Changes
figure_config.plot_retro2( fig1.axes, t_s, y3(:,state), xhat3(:,state), {'g-','linewidth',3.0});
% Add Node
figure_config.plot_retro2( fig1.axes, t_s, y1(:,state), xhat1(:,state), {'b-','linewidth',0.8});
% Add Controller
figure_config.plot_retro2( fig1.axes, t_s, y2(:,state), xhat2(:,state), {'r-','linewidth',0.8});
% Remove Controller
figure_config.plot_retro2( fig1.axes, t_s, y4(:,state), xhat4(:,state), {'k-','linewidth',0.8});
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
%% Calculation vhat
vhat1 = lsim( sys_env_o1, w1, t_s);
vhat2 = lsim( sys_env_o1, w2, t_s);
vhat3 = lsim( sys_env_o1, w3, t_s);
vhat4 = lsim( sys_env_o1, w4, t_s);
%% InterConnection Signal
fig2 = figure_config.set_figure_retro('InterConnection v',[0,0,900,900],{'v','\hat{v}','\tilde{v}'});
% Not Changes
vtilde3 = figure_config.plot_retro1( fig2.axes, t_s, v3, vhat3, {'g-','linewidth',3.0});
% Add Node
vtilde1 = figure_config.plot_retro1( fig2.axes, t_s, v1, vhat1, {'b-','linewidth',0.8});
% Add Controller
vtilde2 = figure_config.plot_retro1( fig2.axes, t_s, v2, vhat2, {'r-','linewidth',0.8});
% Remove Controller
vtilde4 = figure_config.plot_retro1( fig2.axes, t_s, v4, vhat4, {'k-','linewidth',0.8});

ex2_ori = fig2.axes.ax1.Children(end);
ex2_add = fig2.axes.ax1.Children(end-1);
ex3_ori = fig2.axes.ax1.Children(end-2);
ex1_ori = fig2.axes.ax1.Children(end-3);
%ex2_int = fig2.axes.ax1.Children(end-4);
legend(fig2.axes.ax1,[ex2_ori,ex2_add,ex3_ori,ex1_ori],...
    'Original(Extend2)','AddNode(Extend2)','Original(Extend3)','Original(Extend1)',...
    'location','best')
%% Calculation dvhat
start_t = 800;
end_t = 900;
start_t = start_t/Ts_s+1;
end_t = end_t/Ts_s+1;

[dvhat1,sys1] = pem_gradient_continuous(tf(sys_env_o1), w1, t_s);
[dvhat2,sys2] = pem_gradient_continuous(tf(sys_env_o1), w2, t_s);
[dvhat3,sys3] = pem_gradient_continuous(tf(sys_env_o1), w3, t_s);
[dvhat4,sys4] = pem_gradient_continuous(tf(sys_env_o1), w4, t_s);

fi3 = figure('Name','garadient');
plot(mean(dvhat1(start_t:end_t,:)),'b')
hold on; grid on
plot(mean(dvhat2(start_t:end_t,:)),'r')
plot(mean(dvhat3(start_t:end_t,:)),'g')
plot(mean(dvhat4(start_t:end_t,:)),'k')
legend('add node','add controller','not change','remove controller' ,'location','best')

maxlag = 15;
correlation1 = correlation_cal(v1(start_t:end_t), dvhat1(start_t:end_t,:), maxlag);
correlation2 = correlation_cal(v2(start_t:end_t), dvhat2(start_t:end_t,:), maxlag);
correlation3 = correlation_cal(v3(start_t:end_t), dvhat3(start_t:end_t,:), maxlag);
[correlation4,lgs] = correlation_cal(v4(start_t:end_t), dvhat4(start_t:end_t,:), maxlag);
params = tf(sys_env_o1);
params = numel(params.Numerator{:}) + (numel(params.Denominator{:})-1);
% for itr = 1 : params
%     figure('Name',strcat('correlation of parameter',num2str(itr),' between \tilde{v} & \partial \hat{v}'))
%     plot(lgs, correlation3(itr,:))
%     hold on; grid on;
%     plot(lgs, correlation1(itr,:))
%     plot(lgs, correlation2(itr,:))
%     plot(lgs, correlation4(itr,:))
%     legend('original','add node','add controller','remove controller','location','best')
% end
fig4 = figure('Name','correlation between \tilde{v} & \partial \hat{v}');
stem( correlation3(:,maxlag+1),'g')
hold on; grid on;
stem( correlation1(:,maxlag+1),'b')
stem( correlation2(:,maxlag+1),'r')
stem( correlation4(:,maxlag+1),'k')
legend('original','add node','add controller','remove controller','location','best')
%%
% figure('Name','update system')
% bode(sys_env_o1)
% hold on
% bode(sys_env_a)
% bode(sys1)
% bode(sys2)
% bode(sys3)
% bode(sys4)
% legend('original','addnodesystem','add node','add controller','not change','remove controller','location','best')
%% local function
