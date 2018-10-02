clear 
close all
%% Generate Network
seed = 10;
Node_number = 10;
n_ori = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% control & noise Node
c_n = 1;
n_n = [2,3];
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
n_ori.add_controller( 2, Q, R);
sys_ori_c1 = n_ori.get_sys_controlled(sys_ori);
%% Generate v & w
N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;
%% fitting
maxitr = 1;
state = 4;
in = 1;
out = 1;
% get memory
Theta1 = zeros(maxitr, 2*state+1);
Mag = zeros(length(omega), maxitr);
Pha = zeros(length(omega), maxitr);

parfor_progress(maxitr);
%%%%%% init sys (coloum vector)
% rng('shuffle')

% init_params = 1e-3;

for itr = 1:maxitr
    d = randn(N,2+numel(n_n)*2);
    d(:,1:2) = (d(:,1:2)+0.5)*id_in_p;
    d(:,3:end) = d(:,3:end)*noise_p;
    % Response of v&w 
%     rng(28);
    rng('shuffle')
    v = lsim(sys_ori_c1(ob_v_p, cat(2,ID_in_p,Noise)), d, t);
    w = lsim(sys_ori_c1(ob_w_p, cat(2,ID_in_p,Noise)), d, t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % best cost calculate
    Best_cost = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out), sys_local_vw, [0 1]));
    Best_cost.set_sys(balred(sys_env ,state));
    best_cost = Best_cost.eval_func(t, [w, v], zeros(N, 1));
    fprintf('Best Cost is %e.', best_cost);
    fprintf('\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% Choice ID SS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For Identification
%     m = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out), sys_local_vw, [0 1]));
%     m = model_ss(gen_ss_rectifier(gen_ss_canonical(state, in, out), sys_local_vw, [0 1]));
    m = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, in, out), sys_local_vw, [0 1]));
%     m = model_ss(gen_ss_rectifier(gen_ss_canonical_zero(state, in, out), sys_local_vw, [0 1]));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%w%%%%%%%%%%%%%%
    % Open Loop ID
    data = iddata(v,w,Ts);
    sys_armax = armax(data, [state, state+1, state, 0]);
    sys_arx = arx(data, [state, state+1, 0]);
    sys_oe = oe(data, [state+1, state, 0]);
    
%     m.set_sys(balred(sys_env,state));
%     m.set_sys(balred(ss(d2c(sys_arx ,'foh')),state));
%     theta = m.get_params_all();
    
%%%%%%%%%%%%%%% fitting
    m.max_iter = 10000;
    % optimizer
    [theta ,J] = m.fit(t, [w, v], zeros(N, 1));
    model = m.gen_ss.gen_ss.get_sys();
    [num, den] = tfdata(model,'v');
    Theta1(itr, :) = [den(2:end), num];
    % model character
    [mag, pha] = bode(model, omega);
    mag = squeeze(mag(1,1,:));
    pha = squeeze(pha(1,1,:));
    Mag(:, itr) = mag;
    Pha(:, itr) = pha;

    parfor_progress();
end
parfor_progress(0);

%% simlation
Ts_s = 0.01;
Tf = 200;
t_s = 0:Ts_s:Tf;
sim_noise = randn(length(t_s), 2);
%% controller parameter
Q = kron(eye(1*numel(1)),diag([1,1000]));
R = kron(eye(1*numel(1)),diag([1e-3]));
%% add controller (ID) & simlation
n_ori.controllers = {};
n_ori.add_controller( c_n, model, Q, R);
sys_ori_ID = n_ori.get_sys_controlled(sys_ori);
y_ID = lsim(sys_ori_ID(ob_y_p, ID_in_p), sim_noise, t_s);
xhat_ID = lsim(sys_ori_ID(ob_xhat_p, ID_in_p), sim_noise, t_s);
%% add controlller (arx_sys) & simlation
n_ori.controllers = {};
n_ori.add_controller(c_n, ss(d2c(sys_armax,'foh')), Q, R);
sys_ori_armax = n_ori.get_sys_controlled(sys_ori);
y_armax = lsim(sys_ori_armax(ob_y_p, ID_in_p), sim_noise, t_s);
xhat_armax = lsim(sys_ori_armax(ob_xhat_p, ID_in_p), sim_noise, t_s);
%% add controlller (arx_sys) & simlation
n_ori.controllers = {};
n_ori.add_controller(c_n, ss(d2c(sys_arx,'foh')), Q, R);
sys_ori_armax = n_ori.get_sys_controlled(sys_ori);
y_arx = lsim(sys_ori_armax(ob_y_p, ID_in_p), sim_noise, t_s);
xhat_arx = lsim(sys_ori_armax(ob_xhat_p, ID_in_p), sim_noise, t_s);
%% add controlller (arx_sys) & simlation
n_ori.controllers = {};
n_ori.add_controller(c_n, ss(d2c(sys_oe,'foh')), Q, R);
sys_ori_armax = n_ori.get_sys_controlled(sys_ori);
y_oe = lsim(sys_ori_armax(ob_y_p, ID_in_p), sim_noise, t_s);
xhat_oe = lsim(sys_ori_armax(ob_xhat_p, ID_in_p), sim_noise, t_s);
%% add controlller (arx_sys) & simlation
n_ori.controllers = {};
n_ori.add_controller(c_n, balred(sys_env, state), Q, R);
sys_ori_ideal = n_ori.get_sys_controlled(sys_ori);
y_ideal = lsim(sys_ori_ideal(ob_y_p, ID_in_p), sim_noise, t_s);
xhat_ideal = lsim(sys_ori_ideal(ob_xhat_p, ID_in_p), sim_noise, t_s);
%% drawing
H1 = figure_config.set_figure_retro('theta');
% Not tilde, Hat
ytilde_ID_t = figure_config.plot_retro2(H1.axes, t_s, y_ID(:,1), xhat_ID(:,1), {'r-', 'linewidth', 1.0});
ytilde_arx_t = figure_config.plot_retro2(H1.axes, t_s, y_arx(:,1), xhat_arx(:,1), {'b-', 'linewidth', 1.0});
ytilde_armax_t = figure_config.plot_retro2(H1.axes, t_s, y_armax(:,1), xhat_armax(:,1), {'g-', 'linewidth', 1.0});
ytilde_oe_t = figure_config.plot_retro2(H1.axes, t_s, y_oe(:,1), xhat_oe(:,1), {'k-', 'linewidth', 1.0});
ytilde_ideal_t = figure_config.plot_retro2(H1.axes, t_s, y_ideal(:,1), xhat_ideal(:,1), {'m-', 'linewidth', 1.0});

H2 = figure_config.set_figure_retro('omega');
ytilde_ID_o = figure_config.plot_retro2(H2.axes, t_s, y_ID(:,2), xhat_ID(:,2), {'r-', 'linewidth', 1.0});
ytilde_arx_o = figure_config.plot_retro2(H2.axes, t_s, y_arx(:,2), xhat_arx(:,2), {'b-', 'linewidth', 1.0});
ytilde_armax_o = figure_config.plot_retro2(H2.axes, t_s, y_armax(:,2), xhat_armax(:,2), {'g-', 'linewidth', 1.0});
ytilde_oe_o = figure_config.plot_retro2(H2.axes, t_s, y_oe(:,2), xhat_oe(:,2), {'k-', 'linewidth', 1.0});
ytilde_ideal_o = figure_config.plot_retro2(H2.axes, t_s, y_ideal(:,2), xhat_ideal(:,2), {'m-', 'linewidth', 1.0});


figure('Name','Bode')
bode(sys_env,balred(sys_env, state),sys_arx,sys_armax,sys_oe,model)
legend('original','balred','ARX','ARMAX','OE','ID')

