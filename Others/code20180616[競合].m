% Environment Identification
% 2 indetification input
%% initiallize workspace
clear 
close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 4;
seed = 10;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;

%% Set I/O node 
c_n = 1;
n_n = 4;
name = strcat('_Cnode',num2str(c_n),'_Nnode',num2str(n_n));
%% add I/O port for identificaiton
sys_org = n.get_sys();
for i =  1 : Node_number
    sys_org = n.add_io(sys_org,i, strcat('node',num2str(i)));
end

%% Set Identificaiton Input
N = 50000;
Ts = 0.01;
t = (0:N-1)'*Ts;
cn1 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2);
d1 = cn1();
d1 = d1*1;
cn2 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2);
d2 = cn2();
d2 = d2*1;

%% Response simulation
lsim_type = 'foh';
R = 25;
% Model dimension
dim = 6;
opt = oeOptions;
opt.Display = 'on';
%opt.SearchOptions.MaxIterations = (dim*2+1)*100;
%opt.SearchOptions.Tolerance = 1e-4;
% Response of v&w in single id input (local)
v_l = lsim(sys_org({strcat('v_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), d1, t,lsim_type);
w_l = lsim(sys_org({strcat('w_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), d1, t,lsim_type);
data_l_o = iddata(v_l,w_l,Ts);
data_l = resample(data_l_o,1,R);
init_sys_l = oe(data_l,[dim,dim,0],opt);
model_l_p = d2c(init_sys_l);
model_l = ss(model_l_p);
% Response of v&w in single id input (environment)
v_e = lsim(sys_org({strcat('v_node',num2str(c_n))}, {strcat('d_node',num2str(n_n))}), d2, t,lsim_type);
w_e = lsim(sys_org({strcat('w_node',num2str(c_n))}, {strcat('d_node',num2str(n_n))}), d2, t,lsim_type);
data_e_o = iddata(v_e,w_e,Ts);
data_e = resample(data_e_o,1,R);
init_sys_e = oe(data_e,[dim,dim,0],opt);
model_e_p = d2c(init_sys_e);
model_e = ss(model_e_p);
% Response of v&w in double id input (local&environment)
v_le = lsim(sys_org({strcat('v_node',num2str(c_n))}, {strcat('d_node',num2str(c_n)),strcat('d_node',num2str(n_n))}), [d1,d2], t,lsim_type);
w_le = lsim(sys_org({strcat('w_node',num2str(c_n))}, {strcat('d_node',num2str(c_n)),strcat('d_node',num2str(n_n))}), [d1,d2], t,lsim_type);
data_le_o = iddata(v_le,w_le,Ts);
data_le = resample(data_le_o,1,R);
init_sys_le = oe(data_le,[dim,dim,0],opt);
model_le_p = d2c(init_sys_le);
model_le = ss(model_le_p);

%% character area
[sys_local, sys_env] = n.get_sys_local(c_n);
sys_local_vw = ss(sys_local.A,sys_local.B(:,sys_local.InputGroup.v),sys_local.C(sys_local.OutputGroup.w,:),[]);
FigBode=figure('Name',strcat('BodeDiagram',name));
bode(sys_env,'k:')
hold on
bode(model_l,'r')
bode(model_e,'b')
bode(model_le,'g')
legend('original','Local','Env','Local&Env')

%% Add Controller (Extend Retro fit)
Q = diag([1,1000]);
R = diag([1e-3]);
% Obj.controllers initialize
% Single ID input (local)
n.controllers = {};
n.add_controller(c_n,model_l, Q, R);
sys_cL = n.get_sys_controlled(sys_org);
sys_K_cL = n.controllers{1}.sys_K;
% Single ID input (Env)
n.controllers = {};
n.add_controller(c_n,model_e, Q, R);
sys_cE = n.get_sys_controlled(sys_org);
sys_K_cE = n.controllers{1}.sys_K;
% Single ID input (local&Env)
n.controllers = {};
n.add_controller(c_n,model_le, Q, R);
sys_cLE = n.get_sys_controlled(sys_org);
sys_K_cLE = n.controllers{1}.sys_K;

%% response
s_t = 300;%simulation_time
time = 0:0.01:s_t-0.01;
time = time';
% Real Response
y = impulse(sys_org({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), time);
yL = impulse(sys_cL({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), time);
yE = impulse(sys_cE({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), time);
yLE = impulse(sys_cLE({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), time);
% Designed Response
yL_desire = impulse(sys_K_cL('y','x'),time);
yE_desire = impulse(sys_K_cE('y','x'),time);
yLE_desire = impulse(sys_K_cLE('y','x'),time);
% Noise Response
yL_n = lsim(sys_cL({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(n_n))}), d2(1:numel(time),:), time);
yE_n = lsim(sys_cE({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(n_n))}), d2(1:numel(time),:), time);
yLE_n = lsim(sys_cLE({strcat('y_node',num2str(c_n))}, {strcat('d_node',num2str(n_n))}), d2(1:numel(time),:), time);


check = 1;
norm_number = 2;
Fig2 = figure('Name',strcat('original',name));
plot(time,y(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y(:,1,check),norm_number)),num2str(norm(y(:,2,check),norm_number)),'location','best')

Fig3 = figure('Name',strcat('Extend(LocalModeling)_wo_noise',name));
subplot(3,1,1)
plot(time,yL(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yL(:,1,check),norm_number)),num2str(norm(yL(:,2,check),norm_number)),'location','best')
subplot(3,1,2)
plot(time,yL_desire(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yL_desire(:,1,check),norm_number)),num2str(norm(yL_desire(:,2,check),norm_number)),'location','best')
subplot(3,1,3)
re = yL-yL_desire;
plot(time,re(:,:,check),'LineWidth',1.5)
legend(num2str(norm(re(:,1,check),norm_number)),num2str(norm(re(:,2,check),norm_number)),'location','best')

Fig4 = figure('Name',strcat('Extend(LocalModeling)_w_noise',name));
subplot(3,1,1)
plot(time,yL(:,:,check)+yL_n(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yL(:,1,check),norm_number)),num2str(norm(yL(:,2,check),norm_number)),'location','best')
subplot(3,1,2)
plot(time,yL_desire(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yL_desire(:,1,check),norm_number)),num2str(norm(yL_desire(:,2,check),norm_number)),'location','best')
subplot(3,1,3)
re = yL+yL_n-yL_desire;
plot(time,re(:,:,check),'LineWidth',1.5)
legend(num2str(norm(re(:,1,check),norm_number)),num2str(norm(re(:,2,check),norm_number)),'location','best')

Fig5 = figure('Name',strcat('Extend(EnvModeling)_wo_noise',name));
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

Fig6 = figure('Name',strcat('Extend(EnvModeling)_w_noise',name));
subplot(3,1,1)
plot(time,yE(:,:,check)+yE_n(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yE(:,1,check),norm_number)),num2str(norm(yE(:,2,check),norm_number)),'location','best')
subplot(3,1,2)
plot(time,yE_desire(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yE_desire(:,1,check),norm_number)),num2str(norm(yE_desire(:,2,check),norm_number)),'location','best')
subplot(3,1,3)
re = yE+yE_n-yE_desire;
plot(time,re(:,:,check),'LineWidth',1.5)
legend(num2str(norm(re(:,1,check),norm_number)),num2str(norm(re(:,2,check),norm_number)),'location','best')

Fig7 = figure('Name',strcat('Extend(Local&Env)_wo_noise',name));
subplot(3,1,1)
plot(time,yLE(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yLE(:,1,check),norm_number)),num2str(norm(yLE(:,2,check),norm_number)),'location','best')
subplot(3,1,2)
plot(time,yLE_desire(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yLE_desire(:,1,check),norm_number)),num2str(norm(yLE_desire(:,2,check),norm_number)),'location','best')
subplot(3,1,3)
re = yLE-yLE_desire;
plot(time,re(:,:,check),'LineWidth',1.5)
legend(num2str(norm(re(:,1,check),norm_number)),num2str(norm(re(:,2,check),norm_number)),'location','best')

Fig8 = figure('Name',strcat('Extend(Local&Env)_w_noise',name));
subplot(3,1,1)
plot(time,yLE(:,:,check)+yLE_n(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yLE(:,1,check),norm_number)),num2str(norm(yLE(:,2,check),norm_number)),'location','best')
subplot(3,1,2)
plot(time,yLE_desire(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yLE_desire(:,1,check),norm_number)),num2str(norm(yLE_desire(:,2,check),norm_number)),'location','best')
subplot(3,1,3)
re = yLE+yLE_n-yLE_desire;
plot(time,re(:,:,check),'LineWidth',1.5)
legend(num2str(norm(re(:,1,check),norm_number)),num2str(norm(re(:,2,check),norm_number)),'location','best')

%% 
%
saveas(FigBode,FigBode.Name,'fig');
saveas(Fig2,Fig2.Name,'fig');
saveas(Fig3,Fig3.Name,'fig');
saveas(Fig4,Fig4.Name,'fig');
saveas(Fig5,Fig5.Name,'fig');
saveas(Fig6,Fig6.Name,'fig');
saveas(Fig7,Fig7.Name,'fig');
saveas(Fig8,Fig8.Name,'fig');
%}