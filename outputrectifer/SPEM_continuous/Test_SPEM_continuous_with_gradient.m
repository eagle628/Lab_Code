close all
clear

P = tf([9,0],[1,1,9]);
K = tf(1,[1,8]);

Ts = 0.1;

S = 1/(1+P*K);
G_yr = P*S;
G_ur = S;

G_yd = S;
G_ud = -K*S;

N = 100000;
t = 0:Ts:Ts*(N-1);t = t';
r = (rand(N,1)-0)*1;
% d = (rand(N,1)-0.5)*1;
cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',1);
d = cn3();
clear cn3

u = lsim(G_ur,r,t) + lsim(G_ud,d,t);
y = lsim(G_yr,r,t) + lsim(G_yd,d,t);

data = iddata(y,u,Ts);
opt = armaxOptions;
d_ID = armax(data,[2,3,2,0]);
% d_ID = arx(data,[2,3,0]);
% d_ID = arx(data,[2,2,1]);

% d_ID = oe(data,[3,2,0]);



% optimiztion option
opt = optimoptions('lsqnonlin');
opt.UseParallel = true;
opt.Display = 'iter-detailed';
opt.FiniteDifferenceType = 'central';
opt.FunctionTolerance = 1e-12;
opt.OptimalityTolerance = 1e-10;
% opt.StepTolerance = 1e-6;
opt.MaxFunctionEvaluations = 5000;
opt.SpecifyObjectiveGradient = true;

% Assign
dim = 2;
sym_gra = SPEM_continuous_gradient(dim,K);
cost_func  = @(params)SPEM_continuous_cost_func(dim,params,data,K,[],sym_gra);

init_sys = tf(d2c(d_ID));
init_params = [init_sys.Denominator{:}(2:end),init_sys.Numerator{:}];

% optimization
params = lsqnonlin(cost_func,init_params,[],[],opt);
A = params(1:dim); 
B = params(dim+1:end);

final_model = tf(B,[1,A]);

% Bode
bode(P)
hold on
bode(d_ID)
bode(final_model)

legend('Target','Init','SPEM')