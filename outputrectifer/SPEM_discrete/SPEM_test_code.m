% function [final_model,init_sys] = SPEM_test_code( n, c_n, n_n, id_in_p, noise_p, id_method, sim_seed)
clear 
close all
%% Generate Network
seed = 3;
Node_number = 2;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;
n.plot()
%% control & noise Node
c_n = 1;
n_n = [2];
%% signal power
id_in_p = 1;
noise_p = 1;
%% simlationi noise seed
sim_seed = 10;
%% init system method
id_method = 'ARX';
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
% cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2+numel(n_n)*2,'RandomStream','mt19937ar with seed','Seed',sim_seed);
% d = cn3();
rng('shuffle')
d = randn(N,2+numel(n_n)*2);
d(:,1:2) = (d(:,1:2)+0.5)*id_in_p;
d(:,3:end) = d(:,3:end)*noise_p;
clear cn3;
% simulation
% R = 10;
% Response of v&w 
v = lsim(sys_org(ob_v, cat(2,ID_in,Noise)), d, t);
w = lsim(sys_org(ob_w, cat(2,ID_in,Noise)), d, t);
data = iddata(v,w,Ts);
%     data = resample(data,1,1);
%% init system
dim = 2*(Node_number-1);
init_sys = arx(data,[dim,dim+1,0]);
% init_sys = armax(data,[dim,dim+1,dim,0]);
% init_sys = oe(data,[dim+1,dim,0]);
%% Stabilized Prediction Error Method
sys_rect = sys_local({'omega','w'},{'v'}); 
sys_rect = c2d(sys_rect,data.Ts,'foh');


%%%%%%%%%%%%%%%% lsqnonlin
% init_TF = tf(init_sys);
% % init_TF = tf(c2d(sys_env,data.Ts));
% init_params = [init_TF.Denominator{:}(2:end),init_TF.Numerator{:}];
% 
% % optimiztion option
% opt = optimoptions('lsqnonlin');
% opt.UseParallel = true;
% opt.Display = 'iter-detailed';
% opt.FiniteDifferenceType = 'central';
% opt.FunctionTolerance = 1e-7;
% opt.OptimalityTolerance = 1e-7;
% opt.StepTolerance = 1e-7;
% opt.MaxFunctionEvaluations = 5000;
% opt.DiffMaxChange = 1e-10;
% % obj.option.SpecifyObjectiveGradient = true;
% 
% % Assign
% cost_func = @(params)Hanse_cost_func(params,data,-sys_rect({'w'}),dim);
% [final_params,resnorm,residual,exitflag,output] = lsqnonlin(cost_func,init_params,[],[],opt);
% %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% partcle swarm
M = 30;
rng('shuffle')
nvar = dim*2+1;
init_params = zeros(M,nvar);
for itr = 1:M
    stable = 0;
    while ~stable
        rand_params = randn(1,2*dim);
        rand_sys = tf([0,rand_params(dim+1:end)],[1,rand_params(1:dim)],data.Ts);
        rand_sys = rand_sys * tf([1,-1],1,data.Ts);
        all = rand_sys*sys_rect({'w'})/(1-rand_sys*sys_rect({'w'}));
        all = minreal(all,[],false);
        stable = isstable(all);
    end
    [num,den] = tfdata(rand_sys,'v');
    init_params(itr,:) = [den(2:end),num];
end

opt = optimoptions('particleswarm'); 
opt.Display = 'iter';
opt.SwarmSize = M;
opt.MaxStallIterations = 20;
opt.InitialSwarmMatrix = init_params;
opt.UseParallel = true;
cost_func_pso = @(params)sum(sum((Hanse_cost_func(params,data,-sys_rect({'w'}),dim)).^2));
[final_params,fval,exitflag,output] = particleswarm(cost_func_pso, nvar, [], [], opt);
%%%%%%%%%%%%%%%

A = final_params(1:dim); 
B = final_params(dim+1:end);

final_model = tf(B,[1,A],data.Ts);

compare = armax(data,[dim,dim+1,dim,0]);
%%
figure
bode(sys_env,init_sys,final_model)
hold on
init_sys_set = cell(1,M);
for itr = 1:M
den = [1,init_params(itr,1:dim)];
num = init_params(itr,dim+1:end);
sys = tf(num,den,data.Ts);
init_sys_set{itr} = sys;
end
for itr =1:M
bode(init_sys_set{itr},'g:')
end

%%%%%%%%%%%%%%%%%%
% H = figure_config.set_figure_bode('');
% 
% min = -6;
% max =  3;
% [mag_env,pha_env,omega_env] = bode(sys_env,{10^min,10^max});
% mag_env=squeeze(mag_env(1,1,:));
% pha_env=squeeze(pha_env(1,1,:));
% figure_config.plot_bode(H.axes,omega_env,mag_env,pha_env,{'b-','linewidth',3.0})
% 
% [mag_env,pha_env,omega_env] = bode(final_model,{10^min,10^max});
% mag_env=squeeze(mag_env(1,1,:));
% pha_env=squeeze(pha_env(1,1,:));
% figure_config.plot_bode(H.axes,omega_env,mag_env,pha_env,{'r-','linewidth',2.0})
% 
% [mag_env,pha_env,omega_env] = bode(init_sys,{10^min,10^max});
% mag_env=squeeze(mag_env(1,1,:));
% pha_env=squeeze(pha_env(1,1,:));
% figure_config.plot_bode(H.axes,omega_env,mag_env,pha_env,{'g:','linewidth',2.0})

% [mag_env,pha_env,omega_env] = bode(compare,{10^min,10^max});
% mag_env=squeeze(mag_env(1,1,:));
% pha_env=squeeze(pha_env(1,1,:));
% figure_config.plot_bode(H.axes,omega_env,mag_env,pha_env,{'k-','linewidth',2.0})
