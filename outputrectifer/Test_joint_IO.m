%close all
clear
%init_sys_set = cell(1,1000);
%final_model_set = cell(1,1000);
%for itr = 1 : 1000
    
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 30;
seed = 10;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;

%% Set I/O node 
c_n = 1;
n_n = [2,3,4];
%% Signal Power
id_in_p = 1;
noise_p = 0.1;
%% Identification Method
id_method = 'OE';

%% Figure Name
ob_v  = {strcat('v_node',num2str(c_n))};
ob_w  = {strcat('w_node',num2str(c_n))};
ob_y  = {strcat('y_node',num2str(c_n))};
ID_in = {strcat('d_node',num2str(c_n))};
Noise = cell(1,numel(n_n));
for i = 1 : numel(n_n)
    Noise(i) = {strcat('d_node',num2str(n_n(i)))};
end

name_c = '';
for i = 1 : numel(c_n)
    name_c = strcat(name_c,num2str(c_n(i)),'+');
end
name_n = '';
for i = 1 : numel(n_n)
    name_n = strcat(name_n,num2str(n_n(i)),'+');
end
str = replace(num2str(noise_p),'.','%');
name = strcat(id_method,'_N',num2str(Node_number),'_Cnode',name_c,'_Nnode',name_n,'_IDinput',num2str(id_in_p),'_Noise',str);
%% add I/O port for identificaiton
sys_org = n.get_sys();
for i =  1 : n.N
    sys_org = n.add_io(sys_org,i, strcat('node',num2str(i)));
end
min = -10;
max =  10;
[sys_local, sys_env] = n.get_sys_local(c_n);
[mag_env,phase_env,omega] = bode(sys_env,{10^min,10^max});
mag_env=reshape(mag_env(1,1,:),[],1);
phase_env=reshape(phase_env(1,1,:),[],1);
%% Add Controller (Simple Retro fit)
Q = kron(eye(1*numel(c_n)),diag([1,1000]));
R = kron(eye(1*numel(c_n)),diag([1e-3]));
% Ideal 
n.controllers = {};
n.add_controller(c_n,ss([],[],[],0), Q, R);
controlled_sys_I = n.get_sys_controlled(sys_org);
controlled_sys_I_K = n.controllers{1}.sys_K;
%% Generate v & w
N = 100000;
Ts = 0.01;
t = (0:N-1)'*Ts;
rng();
d = rand(N,(1+numel(n_n))*2);
d(:,1:2) = d(:,1:2)*id_in_p;
d(:,3:end) = (d(:,3:end)-0.5)*noise_p;
% simulation
lsim_type = 'foh';
R = 1;
% Response of v&w 
v = lsim(sys_org(ob_v, cat(2,ID_in,Noise)), d, t,lsim_type);
w = lsim(sys_org(ob_w, cat(2,ID_in,Noise)), d, t,lsim_type);
%data = iddata(v,w,Ts);
