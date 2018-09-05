% Environment Identification
%{
Implements Stabilized Prediction Error Method 
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
n_ori.plot()
%% Controlled node number
c_n = 1;
%% Identification SN
id_SN = 0.1;
id_n = [1,3,4];      % identification noise point
%% Number of Iteration
ITR = 1; 
%% Node & Figure Name
% identificaiton noise point
id_noise_p = cell(1,numel(id_n));
for i = 1 : numel(id_n)
    id_noise_p(i) = {strcat('d_node',num2str(id_n(i)))};
end
% interconnection v,w
ob_v_p  = {strcat('v_node',num2str(c_n))};
ob_w_p  = {strcat('w_node',num2str(c_n))};

%% add I/O port for original system
sys_ori = n_ori.get_sys();
add_io = unique([c_n,id_n]);
for i =  1 : numel(add_io)
    sys_ori = n_ori.add_io(sys_ori,add_io(i), strcat('node',num2str(add_io(i))));
end
%% Environment Character (Original)
% Chose Range Fruquency
min = -6;
max =  6;
[sys_local, sys_env] = n_ori.get_sys_local(c_n);
%% generate v & w
% input noise
N = 100000;
Ts = 0.01;
time = 0:Ts:Ts*(N-1);
time = time';
rng('shuffle');
% sim_seed = randi(1000,1,1);
% cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(id_n),'RandomStream','mt19937ar with seed','Seed',sim_seed);
% noise = cn3();
% clear cn3
noise = randn(N,2*numel(id_n));

sim_noise = noise*0.01;
% main
v = lsim(sys_ori(ob_v_p,id_noise_p), sim_noise, time);
w = lsim(sys_ori(ob_w_p,id_noise_p), sim_noise, time);
data = iddata(v,w,Ts);
%% Stabilized Prediction Error Method (continuous)
sys_rect = sys_local({'omega','w'},{'v'});

init_sysd = arx(data,[6,7,0]);
init_sys = tf(d2c(init_sysd));

init_params = [init_sys.Denominator{:}(2:end),init_sys.Numerator{:}];

% optimiztion option
opt = optimoptions('lsqnonlin');
opt.Algorithm = 'levenberg-marquardt';
opt.Display = 'iter-detailed';
opt.FiniteDifferenceType = 'central';
opt.FunctionTolerance = 1e-10;
opt.StepTolerance = 1e-10;
opt.MaxFunctionEvaluations = 5000;
opt.DiffMaxChange = 1e-10;

% Assign
dim = 6;
cost_func = @(params)SPEM_cost_func2(params,data,sys_rect,dim);

[final_params,resnorm,residual,exitflag,output] = lsqnonlin(cost_func,init_params,[],[],opt);
A = final_params(1:dim); 
B = final_params(dim+1:end);

final_model = tf(B,[1,A]);
%%
bode(final_model,sys_env,init_sys)

