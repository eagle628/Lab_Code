close all
clear

P_s = tf([9,0],[1,1,9]);
K_s = tf(1,[1,8]);

Ts = 0.1;
P = c2d(P_s,Ts);
K = c2d(K_s,Ts);

S = 1/(1+P*K);
G_yr = P*S;
G_ur = S;

G_yd = S;
G_ud = -K*S;

N = 100000;
r = (rand(N,1)-0)*1;
% d = (rand(N,1)-0.5)*1;
cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',1);
d = cn3();
clear cn3

u = lsim(G_ur,r) + lsim(G_ud,d);
y = lsim(G_yr,r) + lsim(G_yd,d);

data = iddata(y,u,Ts);
opt = armaxOptions;
d_ID = armax(data,[2,3,2,0]);
% d_ID = arx(data,[2,3,0]);
% d_ID = arx(data,[2,2,1]);

% d_ID = oe(data,[3,2,0]);



% optimiztion option
opt = optimoptions('lsqnonlin');
opt.Display = 'iter-detailed';
opt.FiniteDifferenceType = 'central';
opt.FunctionTolerance = 1e-12;
opt.OptimalityTolerance = 1e-10;
% opt.StepTolerance = 1e-6;
opt.MaxFunctionEvaluations = 100000;

% Assign
dim = 2;
cost_func  = @(params)Hanse_cost_func(params,data,K,dim);

init_sys = tf(d_ID);
init_params = [init_sys.Denominator{:}(2:end),init_sys.Numerator{:}(2:end)];
% init_params = [P.Denominator{:}(2:end),P.Numerator{:}(2:end)];

% optimization
params = lsqnonlin(cost_func,init_params,[],[],opt);
A = params(1:dim); 
B = params(dim+1:end);

final_model = tf(B,[1,A],Ts);

% Bode
bode(P)
hold on
bode(d_ID)
bode(final_model)

legend('Target','Init','SPEM')