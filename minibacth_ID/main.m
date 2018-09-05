clear 
close all
addpath('C:\Users\Naoya Inoue\Documents\GitHub\Lab_Code\sysID_adam_k\systemID')
s = tf('s');
rng(100)
Ts = 0.1;
sys = 1/(s+1);
K = 1/(s+2);

sys_all = 1/(1+sys*K);

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

N = 1000;
d = lsim(H, randn(N,1)*1);
r = randn(N, 1)*1;
t = (0:N-1)'*Ts;

v = lsim(Gv, d, t) + lsim(Gvr, r, t);
w = lsim(Gw, d, t) + lsim(Gwr, r, t);

m = model_ss(gen_ss_rectifier(gen_ss_tridiag(1, 1, 1), K));
m.add_fixed_params('theta_D_1', 0); % model.m
% m.str_display='none';
m.max_iter = 1e-4;
% m.fit(t, [w, v], zeros(N, 1));
m.fit_adam(t, [w, v], zeros(N, 1), 1e-5);
model1 = tf(m.gen_ss.gen_ss.get_sys());
Theta1(itr, :) = [model1.Denominator{1}(2), model1.Numerator{1}(2)];

parfor_progress();
end
parfor_progress(0);


