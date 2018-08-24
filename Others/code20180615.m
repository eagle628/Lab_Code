% Environment Identification
% other point in noise
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
node = [2,27];
c_n = 2;

name = strcat('Cnode',num2str(node(1)),'_Nnode',num2str(node(2)));
%% add I/O port for identificaiton
sys_org = n.get_sys();
for i =  1 : numel(node)
    sys_org = n.add_io(sys_org,node(i), strcat('Id_node',num2str(i)));
end

%% Set Identificaiton Input
N = 50000;
Ts = 0.01;
d = zeros(N,10);
rng(10);
%d(:,[1,2]) = randn(N, 2)*1;
t = (0:N-1)'*Ts;
noise_power = 0.11;
cn = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',10);
d = cn();
d(:,3:end) = d(:,3:end)*noise_power;
%d(:,2)  = zeros(N,1);

%% Response simulation
lsim_type = 'foh';
R = 25;

% Model dimension
dim = 6 ;
system = [dim,dim+1,0];
opt = oeOptions;
opt.Display = 'on';
opt.SearchOptions.MaxIterations = 300;
opt.SearchOptions.Tolerance = 1e-4;
% Response of v&w for d1
v0 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1'}), d(:,1:2), t,lsim_type);
w0 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1'}), d(:,1:2), t,lsim_type);
data0_o = iddata(v0,w0,Ts);
data0 = resample(data0_o,1,R);
%init_sys0 = arx(data0,[6,7,0],opt);
%init_sys0 = armax(data0,[6,7,6,0],opt);
init_sys0 = oe(data0,[6,6,0],opt);
%init_sys0 = bj(data0,[7,6,6,6,0],opt);
model0_p = d2c(init_sys0);
model0 = ss(model0_p);
%model0 = n4sid(data0,6,'N4weight','MOESP','Feedthrough',1,'DisturbanceModel','estimate','Form','canonical','InitialState','zero');
% Response of v&w for d1 & d2
v1 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1','d_Id_node2'}), d(:,1:4), t,lsim_type);
w1 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1','d_Id_node2'}), d(:,1:4), t,lsim_type);
data1_o = iddata(v1,w1,Ts);
data1 = resample(data1_o,1,R);
%init_sys1 = arx(data1,[6,7,0],opt);
%init_sys1 = armax(data1,[6,7,6,0],opt);
init_sys1 = oe(data1,[6,6,0],opt);
%init_sys1 = bj(data1,[7,6,6,6,0],opt);
model1_p = d2c(init_sys1);
model1 = ss(model1_p);
%model1 = n4sid(data1,6,'N4weight','MOESP','Feedthrough',1,'DisturbanceModel','estimate','Form','canonical','InitialState','zero');

%% Chracter plot Area
[~, sys_env] = n.get_sys_local(c_n);
sys_env_low = balred(sys_env, 6);
FigBode = figure('Name',strcat('Bode_',name));
bode(sys_env,'b:');
hold on; grid on;
bode(model0,'g');
bode(model1,'r');
bode(sys_env_low,'k');

legend('original','Ideal','real','balred')
%% Add Controller (Extend Retro fit)
sys_all = n.get_sys();
sys_all = n.add_io(sys_all,c_n, strcat('node',num2str(c_n)));
Q = diag([1,1]);
R = diag([1e-3]);
% Obj.controllers initialize
% w/o noise　(Only node1 injection noise)
n.controllers = {};
n.add_controller(c_n,sys_env_low, Q, R);
sys_cE = n.get_sys_controlled(sys_all);
sys_K_cE = n.controllers{1}.sys_K;
% w/o noise　( 1point injection noise)
n.controllers = {};
n.add_controller(c_n,model0, Q, R);
sys_c0 = n.get_sys_controlled(sys_all);
sys_K_c0 = n.controllers{1}.sys_K;
% w/o noise　( 2point injection noise)
n.controllers = {};
n.add_controller(c_n,model1, Q, R);
sys_c1 = n.get_sys_controlled(sys_all);
sys_K_c1 = n.controllers{1}.sys_K;

%% response
s_t = 300;%simulation_time
time = 0:0.01:s_t-0.01;
time = time';
y = impulse(sys_all({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), time);
yE = impulse(sys_cE({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), time);
y0 = impulse(sys_c0({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), time);
y1 = impulse(sys_c1({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), time);

yE_desire = impulse(sys_K_cE('y','x'),time);
y0_desire = impulse(sys_K_c0('y','x'),time);
y1_desire = impulse(sys_K_c1('y','x'),time);

noise = zeros(length(time),2);
noise(:,1) = randn(length(time),1);


check = 1;
norm_number = 2;
Fig1 = figure('Name',strcat('original_',name));
plot(time,y(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y(:,1,check),norm_number)),num2str(norm(y(:,2,check),norm_number)),'location','best')

Fig2 = figure('Name',strcat('Extend_balred_node_',name));
subplot(3,1,1)
plot(time,yE(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yE(:,1,check),norm_number)),num2str(norm(yE(:,2,check),norm_number)),'location','best')
subplot(3,1,2)
plot(time,yE_desire(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yE_desire(:,1,check),norm_number)),num2str(norm(yE_desire(:,2,check),norm_number)),'location','best')
subplot(3,1,3)
re = yE-yE_desire;
plot(time,re(:,:,check),'LineWidth',1.5)
legend(num2str(norm(re(:,1,check),norm_number)),num2str(norm(re(:,2,check),norm_number)),'location','best')

Fig3 = figure('Name',strcat('Extend_modeling(ideal)_',name));
subplot(3,1,1)
plot(time,y0(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y0(:,1,check),norm_number)),num2str(norm(y0(:,2,check),norm_number)),'location','best')
subplot(3,1,2)
plot(time,y0_desire(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y0_desire(:,1,check),norm_number)),num2str(norm(y0_desire(:,2,check),norm_number)),'location','best')
subplot(3,1,3)
re = y0-y0_desire;
plot(time,re(:,:,check),'LineWidth',1.5)
legend(num2str(norm(re(:,1,check),norm_number)),num2str(norm(re(:,2,check),norm_number)),'location','best')

Fig4 = figure('Name',strcat('Extend_modeling(noise)_',name));
subplot(3,1,1)
plot(time,y1(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y1(:,1,check),norm_number)),num2str(norm(y1(:,2,check),norm_number)),'location','best')
subplot(3,1,2)
plot(time,y1_desire(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y1_desire(:,1,check),norm_number)),num2str(norm(y1_desire(:,2,check),norm_number)),'location','best')
subplot(3,1,3)
re = y1-y1_desire;
plot(time,re(:,:,check),'LineWidth',1.5)
legend(num2str(norm(re(:,1,check),norm_number)),num2str(norm(re(:,2,check),norm_number)),'location','best')

%% 
saveas(FigBode,FigBode.Name,'fig');
saveas(Fig1,Fig1.Name,'fig');
saveas(Fig2,Fig2.Name,'fig');
saveas(Fig3,Fig3.Name,'fig');
saveas(Fig4,Fig4.Name,'fig');
