clear variables
close all
load('kawaguchi_method')
model1_tf = tf(model1);
rng(10);
n = network_swing_simple(30, [1, 2], [2, 10]*1e-2, 1, [1, 5], 0.1, 10);

n.Adj_ref = n.Adj_ref*0;
Q = diag([1, 1000]);
R = 1e-3;

sys_org = n.get_sys();
sys_org = n.add_io(sys_org,{1, 3}, {'n1', 'n3'});
N = 50000;
Ts = 0.01;
d = randn(N, 2);
d2 = randn(N, 2);
t = (0:N-1)'*Ts;
v = lsim(sys_org({'v_n1'}, 'd_n1'), d, t);
w = lsim(sys_org({'w_n1'}, 'd_n1'), d, t);

[~, sys_env] = n.get_sys_local(1);
sys_env_low = balred(sys_env,6);

% id data
data = iddata(v,w,Ts);
%ARX
na =6;nb = 7;nk = 0;
model_arx = arx(data,[na nb nk],'Display','on','Focus','simulation');
model_arx_ss = ss(model_arx);

%ARMAX
na = 6;nb = 7;nc = 6;nk = 0;
model_armax = armax(data,[na nb nc nk],'Display','on','Focus','simulation');
model_armax_ss = ss(model_armax);

opt = polyestOptions;
opt.Display = 'on';
opt.SearchMethod = 'lm';
%opt.SearchMethod = 'lsqnonlin' ;
opt.SearchOptions.MaxIterations = 100;
opt.Focus = 'prediction';
%sys = pem(data,model_armax,opt);
sys = pem(data,sys_env_low,opt);
cov = getcov(sys);
compare(data,sys,sys_env_low,model1)
%sys1 = pem(data,idss(sys_env_low),opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
n_model = 6;
m = model_cont_std(n_model, 0);%何次モデルにするか
m.set_tf(model_armax);% 低次元化パラメータを初期値にセット
m.lsim_type = 'foh';
m.fit_mymqt(t, w, v);% マルカール法によるパラメータ推定
model2 = minreal(ss(m));% 最小実現
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[-600,300,600,600])
bode(model1)
hold on
bode(sys_env)
bode(model_armax)
bode(sys)

legend('kawaguchi','envirionment','init','pem')

%{
%4SID
nx = 6;
model_4sid = n4sid(data,nx,'Display','on');
bode(model_4sid)
model_pem = pem(data,model_4sid);
bode(model_pem)
%}
