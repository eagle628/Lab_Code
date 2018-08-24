clear 
close all
%% Generate Network
seed = 2;
n = network_swing_simple(3, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;
%% control & noise Node
c_n = 1;
n_n = [2,3];
%% signal power
id_in_p = 1;
noise_p = 1;
%% Node Name
ob_v  = {strcat('v_node',num2str(c_n))};
ob_w  = {strcat('w_node',num2str(c_n))};
ID_in = {strcat('d_node',num2str(c_n))};
Noise = cell(1,numel(n_n));
for i = 1 : numel(n_n)
    Noise(i) = {strcat('d_node',num2str(n_n(i)))};
end
%% add I/O port for identificaiton
sys_org = n.get_sys();
for i =  1 : n.N
    sys_org = n.add_io(sys_org,i, strcat('node',num2str(i)));
end
[sys_local, sys_env] = n.get_sys_local(c_n);
%% Generate v & w
N = 100000;
Ts = 0.01;
t = (0:N-1)'*Ts;
rng('shuffle')
cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2+numel(n_n)*2);
d = cn3();
d(:,1:2) = d(:,1:2)*id_in_p;
d(:,3:end) = d(:,3:end)*noise_p;
clear cn3;
% simulatio
% Response of v&w 
v = lsim(sys_org(ob_v, cat(2,ID_in,Noise)), d, t);
w = lsim(sys_org(ob_w, cat(2,ID_in,Noise)), d, t);
data = iddata(v,w,Ts);
%% init system
dim = 4;
% init_sys = arx(data,[dim,dim+1,0]);
init_sys = armax(data,[dim,dim+1,dim,0]);
% init_sys = oe(data,[dim+1,dim,0]);
%% Stabilized Prediction Error Method
sys_rect = sys_local({'omega','w'},{'v'});

init_TF = tf(d2c(init_sys));
init_params = [init_TF.Denominator{:}(2:end),init_TF.Numerator{:}];

% optimiztion option
opt = optimoptions('lsqnonlin');
opt.UseParallel = true;
opt.Display = 'iter-detailed';
opt.FiniteDifferenceType = 'central';
opt.FunctionTolerance = 1e-10;
opt.OptimalityTolerance = 1e-10;
opt.StepTolerance = 1e-10;
opt.MaxFunctionEvaluations = 5000;
opt.DiffMaxChange = 1e-8;
% opt.SpecifyObjectiveGradient = true;


% Assign
sym_gra = SPEM_continuous_gradient(dim,-sys_rect({'w'}));
cost_func  = @(params)SPEM_continuous_cost_func(dim,params,data,-sys_rect({'w'}),-sys_rect({'omega'}),sym_gra);

% optimization
[final_params,resnorm,residual,exitflag,output] = lsqnonlin(cost_func,init_params,[],[],opt);

A = final_params(1:dim); 
B = final_params(dim+1:end);

final_model = tf(B,[1,A]);

 bode(sys_env,init_sys,final_model)
 legend('Traget','Init','SPEM')
