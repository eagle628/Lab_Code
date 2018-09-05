clear 
close all
s = tf('s');
rng(100)
Ts = 0.1;
sys = (s+5)/(s+1)/(s+10);
K = 1/(s+2);

sys_all = 1/(1+sys*K);

H = tf(1, 1, Ts);

Gv = 1/(1+sys*K);
Gw = K/(1+sys*K);

Gwr = 1/(1+sys*K);
Gvr = sys/(1+sys*K);
%% Plant character
min = -6;
max =  6;
[mag_env, pha_env, omega] = bode(sys, {10^min, 10^max});
mag_env = squeeze(mag_env(1,1,:));
pha_env = squeeze(pha_env(1,1,:));
%% minibatch
maxitr = 4;
state = 2;
in = 1;
out = 1;
% get memory
Theta1 = zeros(maxitr, 2*state+1);
Mag = zeros(length(omega), maxitr);
Pha = zeros(length(omega), maxitr);

parfor_progress(maxitr);
parfor itr = 1:maxitr

N = 1000;
d = lsim(H, randn(N,1)*1);
r = randn(N, 1)*1;
t = (0:N-1)'*Ts;

v = lsim(Gv, d, t) + lsim(Gvr, r, t);
w = lsim(Gw, d, t) + lsim(Gwr, r, t);

m = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, in, out), K));
m.add_fixed_params('theta_D_1', 0); % model.m
% m.str_display='none';
m.max_iter = 1e4;
% m.fit(t, [w, v], zeros(N, 1));
m.fit_adam(t, [w, v], zeros(N, 1), 1e-5);
%   eval_func内部の離散化は，あくまで，シミュレーション用であるので，
%   出てくるパラメータ自体は，連続時間のもの
model1 = m.gen_ss.gen_ss.get_sys();
[num, den] = tfdata(model1,'v');
Theta1(itr, :) = [den(2:end), num];
% model character
[mag, pha] = bode(model1, omega);
mag = squeeze(mag(1,1,:));
pha = squeeze(pha(1,1,:));
Mag(:, itr) = mag;
Pha(:, itr) = pha;

parfor_progress();
end
parfor_progress(0);

%% Drawing
H = figure_config.set_figure_bode();
figure_config.plot_bode(H.axes,omega,mag_env,pha_env,{'r:','linewidth',3.0});

figure_config.plot_bode(H.axes,omega,Mag,Pha,{'b-','linewidth',0.8});

