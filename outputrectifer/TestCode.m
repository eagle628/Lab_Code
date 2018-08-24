%close all
clear
%init_sys_set = cell(1,1000);
%final_model_set = cell(1,1000);
%for itr = 1 : 1000
    
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 4;
seed = 10;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;

%% Set I/O node 
c_n = 3;
n_n = [4];
%% Signal Power
id_in_p = 1;
noise_p = 1;
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
n.add_controller( c_n, Q, R);
controlled_sys_I = n.get_sys_controlled(sys_org);
controlled_sys_I_K = n.controllers{1}.sys_K;
%% Generate v & w
N = 100000;
Ts = 0.01;
t = (0:N-1)'*Ts;
sim_seed = 10;
cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',numel(n_n)*2+2,'RandomStream','mt19937ar with seed','Seed',sim_seed);
d = cn3();
% d(:,1:2) = d(:,1:2)*id_in_p;
% d(:,3:end) = (d(:,3:end)-0.5)*noise_p;
% simulation
lsim_type = 'foh';
R = 10;
% Response of v&w 
v = lsim(sys_org(ob_v, cat(2,ID_in,Noise)), d, t,lsim_type);
w = lsim(sys_org(ob_w, cat(2,ID_in,Noise)), d, t,lsim_type);
data = iddata(v,w,Ts);
data = resample(data,1,R);
% init_sys = arx(data,[6,7,0]);
% init_sys = armax(data,[6,7,6,0]);
%init_sys = oe(data,[7,6,0]);
init_sys = balred(sys_env,6);
%% Stabilized Prediction Error Method
sys_rect = sys_local({'omega','w'},{'v'});
sys_rect = c2d(sys_rect,data.Ts);

init_TF = tf(c2d(init_sys,data.Ts));
% init_TF = tf(init_sys);
init_params = [init_TF.Denominator{:}(2:end),init_TF.Numerator{:}];
%init_params = [0,0,0,0,0,0,-6,0,0,0,0,0,0];
% optimiztion option
% opt = optimoptions('fminunc');
opt = optimoptions('lsqnonlin');
opt.Display = 'iter-detailed';
opt.FiniteDifferenceType = 'central';
opt.FunctionTolerance = 1e-8;
opt.StepTolerance = 1e-8;
opt.MaxFunctionEvaluations = 100000;
opt.DiffMaxChange = 1e-6;
%obj.option.Algorithm = 'trust-region';
%obj.option.SpecifyObjectiveGradient = true;

% Assign
dim = 6;
cost_func = @(params)testcode_cost_func(params,data,sys_rect,dim);

% optimization
% params = fminunc(cost_func,init_params,opt);
params = lsqnonlin(cost_func,init_params,[],[],opt);
A = params(1:dim); 
B = params(dim+1:end);

final_model = tf(B,[1,A],data.Ts);
%%
compare_id = armax(data,init_TF);
%%
figure
bode(sys_env)
hold on
%bode(init_sys)
bode(final_model)
bode(init_sys)
bode(compare_id)

%% 
%init_sys_set(itr) = {init_sys};
%final_model_set(itr) = {final_model};
%end