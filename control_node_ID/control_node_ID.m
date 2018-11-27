close all
clear
%% Generate Network
seed = 3;
Node_number = 3;
n_ori = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% control & noise Node
c_n = 1;
n_n = [2:3];
%% signal power
id_in_p = 1;
noise_p = 0;
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

n_ori.add_controller( 1, Q, R);
sys_ori_c1 = n_ori.get_sys_controlled(sys_ori);
%% Generate v & w
N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;

%% fitting
max_iter = 1;
state = 4;
in = 1;
out = 1;

parfor_progress(max_iter);
%%%%%% init sys (coloum vector)
rng('shuffle')

rng(28)
d = randn(N,2+numel(n_n)*2);
d(:,1:2) = (d(:,1:2)+0.5)*id_in_p*100;
d(:,3:end) = d(:,3:end)*noise_p;

v_c = lsim(sys_ori_c1(ob_v_p, cat(2,ID_in_p,Noise)), d, t);
w_c = lsim(sys_ori_c1(ob_w_p, cat(2,ID_in_p,Noise)), d, t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % best cost calculate
m_Best = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out),...
                    n_ori.controllers{1}.sys_K({'w'},{'v'}), [0 1]));
m_Best.set_sys(balred(sys_env ,state));
best_cost_set = m_Best.eval_func(t, [w_c, v_c], zeros(N, 1));
%     fprintf('Best Cost is %e.', best_cost);
%     fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Choice ID SS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Identification
% m = model_ss(gen_ss_rectifier(gen_ss_canonical(state, in, out), sys_local_vw, [0 1]));
m_c = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, in, out),...
                    n_ori.controllers{1}.sys_K({'w'},{'v'}), [0 1]));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%w%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% fitting
m_c.max_iter = 1e4;
% m.str_display = 'none';
% optimizer
[theta, J] = m_c.fit(t, [w_c, v_c], zeros(N, 1));
% [theta_c, J_] = m_c.fit_adam(t, [w_c, v_c], zeros(N, 1));
sys_ID_c = m_c.gen_ss.gen_ss.get_sys();
% model character
% [Mag_ID(:, itr), Pha_ID(:, itr)] = bode_data(sys_ID_set{itr},omega);
% % % % % 
% % % % % %%%% original system
% % % % % rng(28)
% % % % % v_o = lsim(sys_ori(ob_v_p, cat(2,ID_in_p,Noise)), d, t);
% % % % % w_o = lsim(sys_ori(ob_w_p, cat(2,ID_in_p,Noise)), d, t);
% % % % % 
% % % % % m_o = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, in, out), sys_local_vw, [0, 1]));
% % % % % m_o.max_iter = 1e4;
% % % % % m_o.fit(t, [w_o, v_o], zeros(N, 1));
% % % % % % [theta_o, J_o] = m_o.fit_adam(t, [w_c, v_c], zeros(N, 1));
% % % % % sys_ID_o = m_o.gen_ss.gen_ss.get_sys();

%% Extend response simlation
Ts_s = 0.01;
Tf = 200;
t_s = 0:Ts_s:Tf;
rng(28);
sim_N = length(t_s);
sim_noise = randn(sim_N, 2);

[y_o, ytilde_o] = sim_res(n_ori, sys_ori, sys_ID_o, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, Q, R);
[y_c, ytilde_c] = sim_res(n_ori, sys_ori, sys_ID_c, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, Q, R);

%% local function
function [mag, pha, omega] = bode_data(sys,omega)
    if nargin < 2
        [mag, pha, omega] = bode(sys);
    else
        [mag, pha] = bode(sys, omega);
    end
    mag = squeeze(mag(1,1,:));
    pha = squeeze(pha(1,1,:));
end

function [y, ytilde] = sim_res(n_ori, sys_ori, model, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, Q, R)
    n_ori.controllers = {};
    n_ori.add_controller( 1, model, Q, R);
    sys_ori_ID = n_ori.get_sys_controlled(sys_ori);
    y = lsim(sys_ori_ID(ob_y_p, ID_in_p), sim_noise, t_s);
    ytilde = lsim(sys_ori_ID(ob_xhat_p, ID_in_p), sim_noise, t_s);
end
