clear 
close all
%% Generate Network
seed = 3;
Node_number = 3;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;
n.plot()
%% control & noise Node
c_n = 2;
n_n = [1];
%% signal power
id_in_p = 1;
noise_p = 0.1;
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
% clear cn3;
rng('shuffle')
den = rand(N,2+numel(n_n)*2);
den(:,1:2) = den(:,1:2)*id_in_p;
den(:,3:end) = (den(:,3:end)-0.5)*noise_p;
% simulation
% R = 10;
% Response of v&w 
v = lsim(sys_org(ob_v, cat(2,ID_in,Noise)), den, t);
w = lsim(sys_org(ob_w, cat(2,ID_in,Noise)), den, t);
data = iddata(v,w,Ts);
% data = resample(data,1,5);
%% init system
dim = 2*(Node_number-1);
d_ID = arx(data,[dim,dim+1,0]);
% d_ID = armax(data,[dim,dim+1,dim,0]);
% d_ID = oe(data,[dim+1,dim,0]);
%% Stabilized Prediction Error Method
sys_rect = sys_local({'omega','w'},{'v'});
sys_rect = c2d(sys_rect,data.Ts);

opt = optimoptions('lsqnonlin');
opt.UseParallel = true;
opt.Display = 'iter-detailed';
opt.FiniteDifferenceType = 'central';
opt.FunctionTolerance = 1e-12;
opt.OptimalityTolerance = 1e-10;
opt.StepTolerance = 1e-10;
opt.MaxFunctionEvaluations = 10000;
% opt.SpecifyObjectiveGradient = true;
opt.DiffMaxChange = 1e-8;

% Assign
% sym_gra = SPEM_continuous_gradient(dim,K);
cost_func  = @(params)differentiator_discrete_cost_func(...
                            dim,params,data,-sys_rect({'w'}),-sys_rect({'omega'}),[]);

% init_sys = tf(d2c(d_ID));
% init_params = [init_sys.Denominator{:}(2:end),init_sys.Numerator{:}(1:end-1)];

% zpk_ID = zpk(c2d(sys_env,data.Ts));
zpk_ID = zpk(d_ID);
ZERO = zpk_ID.Z{:};
ZERO_re = real(ZERO);
diff = (1-ZERO_re).^2;
[~,I] = min(diff);
idx = 1:dim;
nidx = setdiff(idx,I);
init_sys_differentiator = zpk(zpk_ID.Z{:}(nidx),zpk_ID.P,zpk_ID.K,data.Ts);
[num,den] = tfdata(init_sys_differentiator,'v');
init_params = [den(2:end),num(2:end)];
init_sys_differentiator = init_sys_differentiator*tf([1,-1],1,data.Ts);

% optimization
params = lsqnonlin(cost_func,init_params,[],[],opt);
A = params(1:dim); 
B = params(dim+1:end);
B = [0,B];

final_model = tf(B,[1,A],data.Ts)*tf([1,-1],1,data.Ts);

figure
bode(sys_env,d_ID,final_model)
legend('original','ARX','SPEM')
%}