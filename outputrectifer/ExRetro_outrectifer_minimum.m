%% initiallize workspace
clear 
close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 10;
seed = 10;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;

%% Set I/O node 
c_n = 1;
n_n = 2;
%% Signal Power
id_in_p = 1;
noise_p = 0.01;
%% Identification Method
id_method = 'ARMAX';
%id_method = 'OE';

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
for i =  1 : numel(c_n)
    sys_org = n.add_io(sys_org,c_n(i), strcat('node',num2str(c_n(i))));
end
for i =  1 : numel(n_n)
    sys_org = n.add_io(sys_org,n_n(i), strcat('node',num2str(n_n(i))));
end

min = -10;
max =  10;
[sys_local, sys_env] = n.get_sys_local(c_n);
[mag_env,phase_env,omega] = bode(sys_env,{10^min,10^max});
mag_env=reshape(mag_env(1,1,:),[],1);
phase_env=reshape(phase_env(1,1,:),[],1);

%% Add Controller ( Retro fit)
Q = kron(eye(1*numel(c_n)),diag([1,1000]));
R = kron(eye(1*numel(c_n)),diag([1e-3]));
%  
n.controllers = {};
n.add_controller(c_n, ss([],[],[],0), Q, R);
controlled_sys = n.get_sys_controlled(sys_org);
controlled_sys_K = n.controllers{1}.sys_K;
%% Rectifer

% Model dimension
dim = 6;
% ID input
N = 100000;
Ts = 0.01;
t = (0:N-1)'*Ts;
% simulation
lsim_type = 'foh';
R = 1;

% SEED
ITR = 1;
rng(6)
rand_s = randi(1000,2,ITR);
i =1;

cn1 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(c_n),'RandomStream','mt19937ar with seed','Seed',rand_s(1,i));
d1 = cn1();
d1 = d1*id_in_p;
cn2 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(n_n),'RandomStream','mt19937ar with seed','Seed',rand_s(2,i));
d2 = cn2();
d2 = d2*noise_p;

% Response of v&w 
v = lsim(sys_org(ob_v, cat(2,ID_in,Noise)), [d1,d2], t,lsim_type);
w = lsim(sys_org(ob_w, cat(2,ID_in,Noise)), [d1,d2], t,lsim_type);
%{
model_arx = arx(iddata(v,w,Ts),[6,7,0]);
den = model_arx.A;
num = model_arx.B;

% Rectifer
G_wv = sys_local({'w'},{'v'});
G_wv = c2d(G_wv,Ts);
G_yv = sys_local({'y'},{'v'});
G_yv = c2d(G_yv,Ts);

% optimiztion option
opt = optimoptions('fminunc');
opt.Display = 'iter-detailed';
opt.FiniteDifferenceType = 'forward';
opt.FunctionTolerance = 1e-6;
opt.StepTolerance = 1e-6;
opt.MaxFunctionEvaluations = 500*(numel(num)+numel(den));

% Assign
params_ini = [den(2:end),num];
cost_func = @(params)out_rect_func(params,w,v,G_wv,G_yv,Ts,6);

% optimization
params = fminunc(cost_func,params_ini,opt);

%
dim = 6;
den = [1,params(1:dim)];
num = params(dim+1:end);
Ge_hat = tf(num,den,Ts);

% Bode
bode(sys_env)
hold on
bode(Ge_hat)
bode(model_arx)
%}