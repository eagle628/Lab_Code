clear 
close all
s = tf('s');
rng(60)
Ts = 0.1;

state = 1;

sys = 1/(s+1);
K = 1/(s+2);

sys_all = 1/(1+sys*K);
pole(sys_all)

H = tf(1, 1, Ts);

Gv = 1/(1+sys*K);
Gw = K/(1+sys*K);

Gwr = 1/(1+sys*K);
Gvr = sys/(1+sys*K);

maxitr = 1;
Theta1 = zeros(maxitr, 2);
Theta2 = zeros(maxitr, 2);
parfor_progress(maxitr);
for itr = 1:maxitr

rng('shuffle')
N = 1000;
d = randn(N,1)*0;
r = randn(N, 1)*1;
t = (0:N-1)'*Ts;

v = lsim(Gv, d, t) + lsim(Gvr, r, t);
w = lsim(Gw, d, t) + lsim(Gwr, r, t);



m = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, 1, 1), K));
% m.str_display='none';
% m.fit(t, [w, v], zeros(N, 1)); 
% m.set_sys(ss(sys))
% J = m.eval_func(t, [w, v], zeros(N, 1), m.gen_ss.gen_ss.theta);
% fprintf('best cost: %e.\n',J);
%%%% init
% init_sys = d2c(arx(iddata(v,w,Ts), [state, state+1, 0]));
% init_params = sys2params_tri(init_sys, state, 1, 1);
%%%% fitting
% m.add_fixed_params('theta_D_1', 0);
theta = 1e-5;
m.max_iter = 5000;
for itr1 = 1 : 10
%     m.fit_adam(t, [w, v], zeros(N, 1), 1e-5);
%     [theta,J_his] = m.fit_Santa_E(t, [w, v], zeros(N, 1), theta, 1e-5, [], [], 1000, @(t)t, 0.5);
    [theta,J_his] = m.fit_Santa_S(t, [w, v], zeros(N, 1), theta, 1e-3, [], 1e-1, 500, @(t)t^2, 0.5);
end

model1 = tf(m.gen_ss.gen_ss.get_sys());
Theta1(itr, :) = [model1.Denominator{1}(2), model1.Numerator{1}(2)];

parfor_progress();
end
parfor_progress(0);

%% 
figure
bode(sys, model1)

%% local function
function init_params = sys2params_tri(init_sys, state, in, out)
    init_sys = canon(balred(ss(init_sys),state),'modal');
    init_params = zeros(3*state-2+state*(in+out)+in*out, 1);
    init_params(1:3*state-2) = [diag(init_sys.A)', diag(init_sys.A, -1)', diag(init_sys.A, 1)'];
    init_params(3*state-2+1:3*state-2+state*in) = reshape(init_sys.B, state*in, 1);
    init_params(3*state-2+state*in+1:3*state-2+state*in+out*state) = reshape(init_sys.C, out*state, 1);
    init_params(3*state-2+state*in+out*state+1:3*state-2+state*in+out*state+in*out) = reshape(init_sys.D, out*in, 1);
end
