clear 
close all
s = tf('s');
rng(100)
Ts = 0.1;
sys = tf([-1.0860,-0.0653,-2.6176,-0.0829,-0.0000],...
            [1.0000,0.0601,3.0582,0.0893,0.9075]);
K = tf([0,0,0.6448],[1,0.0494,0]);
% sys = (s+4)/(s+1)/(s+1);
% K = 1/(s+2);

sys_all = 1/(1+sys*K);

H = tf(1, 1, Ts);

Gvd = 1/(1-sys*K);
Gwd = K/(1-sys*K);

Gwr = 1/(1-sys*K);
Gvr = sys/(1-sys*K);
%% Plant character
min = -6;
max =  6;
[mag_env, pha_env, omega] = bode(sys, {10^min, 10^max});
mag_env = squeeze(mag_env(1,1,:));
pha_env = squeeze(pha_env(1,1,:));
%% minibatch
maxitr = 1;
state = 4;
in = 1;
out = 1;
% get memory
Theta1 = zeros(maxitr, 2*state+1);
Mag = zeros(length(omega), maxitr);
Pha = zeros(length(omega), maxitr);

parfor_progress(maxitr);
for itr = 1:maxitr
rng('shuffle')
N = 1000;
d = lsim(H, randn(N,1)*1);
r = randn(N, 1)*1;
t = (0:N-1)'*Ts;

v = lsim(Gvd, d, t) + lsim(Gvr, r, t);
w = lsim(Gwd, d, t) + lsim(Gwr, r, t);

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

