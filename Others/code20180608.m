% Environment Identification
% other point in noise

% LSで初期値を求め，その後，PEMをfminuncによって求める．
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
N = 50000;
Ts = 0.25;
%d = randn(N, 10);
t = (0:N-1)'*Ts;
noise_power = 1;
cn = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',10);
d= cn()*noise_power;

%% Response simulation

lsim_type = 'foh';

% Model dimension
system = [6,7,0];
% Response of v&w for d1
v0 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1'}), d(:,1:2), t,lsim_type);
w0 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1'}), d(:,1:2), t,lsim_type);
model0 = my_pem_opt(w0,v0,system,t);
%model0 = ss(pem(iddata(v0,w0,Ts),[6,6,0,0,0,0],'Focus','Simulation'));
% Response of v&w for d1 & d2
v1 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1','d_Id_node2'}), d(:,1:4), t,lsim_type);
w1 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1','d_Id_node2'}), d(:,1:4), t,lsim_type);
model1 = my_pem_opt(w1,v1,system,t);
%model1 = ss(pem(iddata(v1,w1,Ts),[6,6,0,0,0,0],'Focus','Simulation'));
% Response of v&w for d1 & d2
v2 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1','d_Id_node2','d_Id_node3'}), d(:,1:6), t,lsim_type);
w2 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1','d_Id_node2','d_Id_node3'}), d(:,1:6), t,lsim_type);
model2 = my_pem_opt(w2,v2,system,t);
%model2 = ss(pem(iddata(v2,w2,Ts),[6,6,0,0,0,0],'Focus','Simulation'));

%% Chracter plot Area
[~, sys_env] = n.get_sys_local(1);
sys_env_low = balred(sys_env, 6);
figure
bode(sys_env);
hold on; grid on;
bode(model0);
bode(model1);
bode(model2);

%% Add Controller (Extend Retro fit)
c_n = 1;
sys_all = n.get_sys();
sys_all = n.add_io(sys_all,c_n, strcat('node',num2str(c_n)));
Q = diag([1,1000]);
R = diag([1e-3]);
% Obj.controllers initialize
% w/o noise　(Only node1 injection noise)
n.controllers = {};
n.add_controller(c_n,sys_env_low, Q, R);
sys_cE = n.get_sys_controlled(sys_all);
% w/o noise　(Only node1 injection noise)
n.controllers = {};
n.add_controller(c_n,model0, Q, R);
sys_c0 = n.get_sys_controlled(sys_all);
% w/o noise　( node1 & node2 injection noise)
n.controllers = {};
n.add_controller(c_n,model1, Q, R);
sys_c1 = n.get_sys_controlled(sys_all);
% w/o noise　( node1 & node2 & node3 injection noise)
n.controllers = {};
n.add_controller(c_n,model2, Q, R);
sys_c2 = n.get_sys_controlled(sys_all);


%% response
s_t = 300;%simulation_time
time = 0:0.01:s_t-0.01;
time = time';
y = impulse(sys_all({'y_node1'}, {'d_node1'}), time);
yE = impulse(sys_cE({'y_node1'}, {'d_node1'}), time);
y0 = impulse(sys_c0({'y_node1'}, {'d_node1'}), time);
y1 = impulse(sys_c1({'y_node1'}, {'d_node1'}), time);
y2 = impulse(sys_c2({'y_node1'}, {'d_node1'}), time);

noise = zeros(length(time),2);
noise(:,1) = randn(length(time),1);


a = 5;
check = 1;
norm_number = 2;
figure
subplot(a,1,1)
plot(time,y(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y(:,1,check),norm_number)),num2str(norm(y(:,2,check),norm_number)),'location','best')
subplot(a,1,2)
plot(time,yE(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yE(:,1,check),norm_number)),num2str(norm(yE(:,2,check),norm_number)),'location','best')
subplot(a,1,3)
plot(time,y0(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y0(:,1,check),norm_number)),num2str(norm(y0(:,2,check),norm_number)),'location','best')
subplot(a,1,4)
plot(time,y1(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y1(:,1,check),norm_number)),num2str(norm(y1(:,2,check),norm_number)),'location','best')
subplot(a,1,5)
plot(time,y2(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y2(:,1,check),norm_number)),num2str(norm(y2(:,2,check),norm_number)),'location','best')
