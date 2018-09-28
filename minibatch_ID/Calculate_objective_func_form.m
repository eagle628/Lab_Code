clear 
close all
%% Generate Network
seed = 3;
Node_number = 30;
n_ori = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% control & noise Node
c_n = 1;
n_n = [2,3];
%% signal power
id_in_p = 1;
noise_p = 0.01;
%% Node Name
ob_v_p  = {strcat('v_node',num2str(c_n))};
ob_w_p  = {strcat('w_node',num2str(c_n))};
ob_y_p  = {strcat('y_node',num2str(c_n))};
ID_in_p = {strcat('d_node',num2str(c_n))};
ob_xhat_p = {'xhat_controlled1'}; % The last number should be Controllers Group number.
Noise = cell(1,numel(n_n));
for i = 1 : numel(n_n)
Noise(i) = {strcat('d_node',num2str(n_n(i)))};
end
%% add I/O port for identificaiton
sys_ori = n_ori.get_sys();
for i =  1 : n_ori.N
sys_ori = n_ori.add_io(sys_ori,i, strcat('node',num2str(i)));
end
[sys_local, sys_env] = n_ori.get_sys_local(c_n);
sys_local_vw = sys_local({'w'},{'v'});
%% Plant character
min = -6;
max =  2;
[mag_env, pha_env, omega] = bode(sys_env, {10^min, 10^max});
mag_env = squeeze(mag_env(1,1,:));
pha_env = squeeze(pha_env(1,1,:));
%% add controller
Q = kron(eye(1*numel(1)),diag([1,1000]));
R = kron(eye(1*numel(1)),diag([1e-3]));
n_ori.add_controller( 2, Q, R);
sys_ori_c1 = n_ori.get_sys_controlled(sys_ori);
%% Generate v & w
N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;
%% minibatch
maxitr = 1;
state = 4;%2*(Node_number-1);
in = 1;
out = 1;
% get memory
Theta1 = zeros(maxitr, 2*state+1);
Mag = zeros(length(omega), maxitr);
Pha = zeros(length(omega), maxitr);

parfor_progress(maxitr);
%%%%%% init sys (coloum vector)
d = randn(N,2+numel(n_n)*2);
d(:,1:2) = (d(:,1:2)+0.5)*id_in_p;
d(:,3:end) = d(:,3:end)*noise_p;
% Response of v&w 
%     rng(28);
rng('shuffle')
v = lsim(sys_ori_c1(ob_v_p, cat(2,ID_in_p,Noise)), d, t);
w = lsim(sys_ori_c1(ob_w_p, cat(2,ID_in_p,Noise)), d, t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% best cost calculate
% Best_cost = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out), sys_local_vw));
% Best_cost.gen_ss.gen_ss.set_sys(sys_env);
% best_cost = Best_cost.eval_func(t, [w, v], zeros(N, 1));
% fprintf('Best Cost is %e.', best_cost);
% fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For Identification
% m = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, in, out), sys_local_vw, [0 1]));
% m = model_ss(gen_ss_rectifier(gen_ss_canonical(state, in, out), sys_local_vw, [0 1]));
m = model_ss(gen_ss_rectifier(gen_ss_canonical_zero(state, in, out), sys_local_vw, [0 1]));
% m = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out), sys_local_vw, [0 1]));
% init
% init_sys = sys_env;
sys_arx = balred(ss(d2c(arx(iddata(v,w,Ts), [state, state+1, 0]), 'foh')), state);
sys_armax = balred(ss(d2c(armax(iddata(v,w,Ts), [state, state+1, state, 0]), 'foh')), state);
sys_oe = balred(ss(d2c(oe(iddata(v,w,Ts), [state+1, state, 0]), 'foh')), state);

% true params
m.gen_ss.gen_ss.set_sys(balred(sys_env,4));
true_params = m.get_params_all();
% Id params
m.gen_ss.gen_ss.set_sys(sys_arx);
arx_params = m.get_params_all();
m.gen_ss.gen_ss.set_sys(sys_armax);
armax_params = m.get_params_all();
m.gen_ss.gen_ss.set_sys(sys_oe);
oe_params = m.get_params_all();

%% estimate objective funciton form
load('sys.mat')
m.gen_ss.gen_ss.set_sys(model);
target_params = m.get_params_all();

init_params = target_params;
% tridiag
number = 1000;
min = -1;
max = 2;
parfor_progress(number);

cost = zeros(1,number);
locale = linspace(min, max, number);
parfor itr = 1 : number
    alpha = locale(itr);
    cost(itr) = m.eval_func(t, [w, v], zeros(N, 1), init_params+alpha.*(true_params-init_params));
    parfor_progress();
end
parfor_progress(0);
id_cost = m.eval_func(t, [w, v], zeros(N, 1), init_params)
best_cost = m.eval_func(t, [w, v], zeros(N, 1), true_params)
%%
close all
figure('Name','ObjectiveFunctionForm')
semilogy(locale,cost)
figure('Name','Bode')
bode(sys_env,sys_arx,sys_armax,sys_oe);
legend('true','ARX','ARMAX','OE')

disp('complete')

%% save
function SAVE(location)
%     location = 'C:\Users\NaoyaInoue\Desktop\figure_set\minibatch\cost_function_form\dim1\01';
    for itr = 1:2
        figure(itr);f = gcf;
        savefig(f, strcat(location,'\',f.Name), 'compact');
    end
    save(strcat(location,'\','data'),...
            'cost','locale',...
            'sys_env','sys_arx','sys_armax','sys_oe',...
            'v','w',...
            'best_cost','id_cost')
end
