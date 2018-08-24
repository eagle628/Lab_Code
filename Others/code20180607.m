% Environment Identification
% other point in noise

% PEMでとりあえずいいや
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
[~, sys_env] = n.get_sys_local(1);
sys_initial = balred(sys_env, 6);

lsim_type = 'foh';

% Model dimension
dim = 6;
na = dim;nb = dim+1;nc = dim;nd = dim;nf = dim;nk = 0;
nn = [na,nb,nc,nd,nf,nk];
% Response of v&w for d1_node
v0 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1'}), d(:,1:2), t,lsim_type);
w0 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1'}), d(:,1:2), t,lsim_type);
data0 = iddata(v0,w0,Ts);
model0 = pem(data0,nn,'Focus','Simulation');
model0_ss = minreal(ss(model0));
% Response of v&w for d1_node&d2_node
v1 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1','d_Id_node2'}), d(:,1:4), t,lsim_type);
w1 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1','d_Id_node2'}), d(:,1:4), t,lsim_type);
data1 = iddata(v1,w1,Ts);
model1 = pem(data1,nn,'Focus','Simulation');
model1_ss = minreal(ss(model1));
% Response of v&w for d1_node&d2_node&d3_node
v2 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1','d_Id_node2','d_Id_node3'}), d(:,1:6), t,lsim_type);
w2 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1','d_Id_node2','d_Id_node3'}), d(:,1:6), t,lsim_type);
data2 = iddata(v2,w2,Ts);
model2 = pem(data2,nn,'Focus','Simulation');
model2_ss = minreal(ss(model2));
% Response of v&w for d1_node&d2_node&d3_node
v3 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1','d_Id_node2','d_Id_node3','d_Id_node4'}), d(:,1:8), t,lsim_type);
w3 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1','d_Id_node2','d_Id_node3','d_Id_node4'}), d(:,1:8), t,lsim_type);
data3 = iddata(v3,w3,Ts);
model3 = pem(data3,nn,'Focus','Simulation');
model3_ss = minreal(ss(model3));



%% Chracter plot Area
figure
bode(sys_env);
hold on; grid on;
bode(model0);
bode(model1);
bode(model2);
bode(model3);

%% Add Controller (Extend Retro fit)
c_n = 1;
sys_all = n.get_sys();
sys_all = n.add_io(sys_all,c_n, strcat('node',num2str(c_n)));
Q = diag([1,1000]);
R = diag([1e-3]);
% Obj.controllers initialize
% w/o noise　(Only node1 injection noise)
n.controllers = {};
n.add_controller(c_n,model0_ss, Q, R);
sys_c0 = n.get_sys_controlled(sys_all);
% w/o noise　( node1,2 injection noise)
n.controllers = {};
n.add_controller(c_n,model1_ss, Q, R);
sys_c1 = n.get_sys_controlled(sys_all);
% w/o noise　( node1,2 injection noise)
n.controllers = {};
n.add_controller(c_n,model2_ss, Q, R);
sys_c2 = n.get_sys_controlled(sys_all);
% w/o noise　( node1,2 injection noise)
n.controllers = {};
n.add_controller(c_n,model3_ss, Q, R);
sys_c3 = n.get_sys_controlled(sys_all);


%% response
s_t = 300;%simulation_time
time = 0:0.01:s_t-0.01;
time = time';
y = impulse(sys_all({'y_node1'}, {'d_node1'}), time);
y0 = impulse(sys_c0({'y_node1'}, {'d_node1'}), time);
y1 = impulse(sys_c1({'y_node1'}, {'d_node1'}), time);
y2 = impulse(sys_c2({'y_node1'}, {'d_node1'}), time);
y3 = impulse(sys_c3({'y_node1'}, {'d_node1'}), time);

noise = zeros(length(time),2);
noise(:,1) = randn(length(time),1);


a = 5;
check = 2;
norm_number = 2;
figure
subplot(a,1,1)
plot(time,y(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y(:,1,check),norm_number)),num2str(norm(y(:,2,check),norm_number)),'location','best')
subplot(a,1,2)
plot(time,y0(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y0(:,1,check),norm_number)),num2str(norm(y0(:,2,check),norm_number)),'location','best')
subplot(a,1,3)
plot(time,y1(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y1(:,1,check),norm_number)),num2str(norm(y1(:,2,check),norm_number)),'location','best')
subplot(a,1,4)
plot(time,y2(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y2(:,1,check),norm_number)),num2str(norm(y2(:,2,check),norm_number)),'location','best')
subplot(a,1,5)
plot(time,y3(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y3(:,1,check),norm_number)),num2str(norm(y3(:,2,check),norm_number)),'location','best')

