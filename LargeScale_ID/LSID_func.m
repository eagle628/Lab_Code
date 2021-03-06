%% NOTE
% ARX, ARMAX of System Idnetification Toolbox cannnot use in PARFOR LOOP.
% So, the part do local function
%% 
% function LSID_func(Node_number, noise_p, n_n, root, max_iter, controller_number, c_n)
% clear 
close all
Node_number=3; noise_p = 0; n_n=[2,3]; root=''; max_iter=1; controller_number=0; c_n=1;
%% Generate Network
seed = 3;
% Node_number = 30;
n_ori = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% control & noise Node
% c_n = 1;
% n_n = [2:3];
%% signal power
id_in_p = 1;
% noise_p = 0;
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

% controller_number = 1;
for itr = 1:controller_number % node1 is target local system
    n_ori.add_controller( itr+1, Q, R);
end
sys_ori_c1 = n_ori.get_sys_controlled(sys_ori);
%% Generate v & w
N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;
%% fitting
% max_iter = 20;
state = 6;
in = 1;
out = 1;

parfor_progress(max_iter);
%%%%%% init sys (coloum vector)
rng('shuffle')

% get memory
Mag_ID = zeros(length(omega), max_iter);
Pha_ID = zeros(length(omega), max_iter);
Mag_arx = zeros(length(omega), max_iter);
Pha_arx = zeros(length(omega), max_iter);
Mag_armax = zeros(length(omega), max_iter);
Pha_armax = zeros(length(omega), max_iter);
Mag_oe = zeros(length(omega), max_iter);
Pha_oe = zeros(length(omega), max_iter);

sys_arx_set = cell(max_iter, 1);
sys_armax_set = cell(max_iter, 1);
sys_oe_set = cell(max_iter, 1);
sys_ID_set = cell(max_iter, 1);

best_cost_set = zeros(max_iter, 1);
ID_cost_set = zeros(max_iter, 1);

parfor itr = 1:max_iter
    % Response of v&w 
%     rng(28);
    d = randn(N,2+numel(n_n)*2);
    d(:,1:2) = (d(:,1:2)+0.5)*id_in_p;
    d(:,3:end) = d(:,3:end)*noise_p;
  
    v = lsim(sys_ori_c1(ob_v_p, cat(2,ID_in_p,Noise)), d, t);
    w = lsim(sys_ori_c1(ob_w_p, cat(2,ID_in_p,Noise)), d, t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % best cost calculate
    m_Best = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out), sys_local_vw, [0 1]));
    m_Best.set_sys(balred(sys_env ,state));
    best_cost_set(itr) = m_Best.eval_func(t, [w, v], zeros(N, 1));
%     fprintf('Best Cost is %e.', best_cost);
%     fprintf('\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% Choice ID SS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For Identification
%     m = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out), sys_local_vw, [0 1]));
%     m = model_ss(gen_ss_rectifier(gen_ss_canonical(state, in, out), sys_local_vw, [0 1]));
    m = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, in, out), sys_local_vw, [0 1]));
%     m = model_ss(gen_ss_rectifier(gen_ss_canonical_zero(state, in, out), sys_local_vw, [0 1]));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%w%%%%%%%%%%%%%%
    % Open Loop ID
    data = iddata(v, w, Ts);
    [sys1, sys2, sys3] = parfor_func(data,state);
    sys_armax_set{itr} = sys1;
    [Mag_armax(:, itr), Pha_armax(:, itr)] = bode_data(sys1,omega);
    sys_arx_set{itr} = sys2;
    [Mag_arx(:, itr), Pha_arx(:, itr)] = bode_data(sys2,omega);
    sys_oe_set{itr} = sys3;
    [Mag_oe(:, itr), Pha_oe(:, itr)] = bode_data(sys3,omega);
%     
%     m.set_sys(balred(sys_env,state));
%     m.set_sys(balred(ss(d2c(sys_arx ,'foh')),state));
%     theta = m.get_params_all();
    
%%%%%%%%%%%%%%% fitting
    m.max_iter = 1e4;
%     m.str_display = 'none';
    % optimizer
    [theta, J] = m.fit(t, [w, v], zeros(N, 1));
    sys_ID_set{itr} = m.gen_ss.gen_ss.get_sys();
    ID_cost_set(itr) = m.eval_func(t, [w, v], zeros(N, 1));
    % model character
%     [mag, pha] = bode(model, omega);
%     mag = squeeze(mag(1,1,:));
%     pha = squeeze(pha(1,1,:));
%     Mag_ID(:, itr) = mag;
%     Pha_ID(:, itr) = pha;
    [Mag_ID(:, itr), Pha_ID(:, itr)] = bode_data(sys_ID_set{itr},omega);
    
    m1 = model_ss(gen_ss_tridiag(state, in, out));
    m1.max_iter = 1e4;
    m1.fit(t, [w, v], zeros(N, 1));
    sys_oe_model = m1.gen_ss.get_sys();
%%%%%%%% for loop end
    parfor_progress();
end
parfor_progress(0);

%% simlation
%%%% Extend response simlation
Ts_s = 0.01;
Tf = 200;
t_s = 0:Ts_s:Tf;
rng(28);
sim_N = length(t_s);
sim_noise = randn(sim_N, 2);
% get memory
y_ID_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
y_arx_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
y_armax_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
y_oe_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
y_ideal_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_ID_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_arx_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_armax_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_oe_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_ideal_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
%%%%%%%% controller parameter
Q = kron(eye(1*numel(1)),diag([1,1000]));
R = kron(eye(1*numel(1)),diag([1e-3]));
for itr = 1 : max_iter
    %%%%%%%% add controller (ID) & simlation
    n_ori.controllers = {};
    n_ori.add_controller( c_n, sys_ID_set{itr}, Q, R);
    sys_ori_ID = n_ori.get_sys_controlled(sys_ori);
    y_ID_set(:, :, itr) = lsim(sys_ori_ID(ob_y_p, ID_in_p), sim_noise, t_s);
    xhat_ID_set(:, :, itr) = lsim(sys_ori_ID(ob_xhat_p, ID_in_p), sim_noise, t_s);
    %%%%%%%% add controlller (armax_sys) & simlation
    n_ori.controllers = {};
    n_ori.add_controller(c_n, ss(d2c(sys_armax_set{itr},'foh')), Q, R);
    sys_ori_armax = n_ori.get_sys_controlled(sys_ori);
    y_armax_set(:, :, itr) = lsim(sys_ori_armax(ob_y_p, ID_in_p), sim_noise, t_s);
    xhat_armax_set(:, :, itr) = lsim(sys_ori_armax(ob_xhat_p, ID_in_p), sim_noise, t_s);
    %%%%%%%% add controlller (arx_sys) & simlation
    n_ori.controllers = {};
    n_ori.add_controller(c_n, ss(d2c(sys_arx_set{itr},'foh')), Q, R);
    sys_ori_arx = n_ori.get_sys_controlled(sys_ori);
    y_arx_set(:, :, itr) = lsim(sys_ori_arx(ob_y_p, ID_in_p), sim_noise, t_s);
    xhat_arx_set(:, :, itr) = lsim(sys_ori_arx(ob_xhat_p, ID_in_p), sim_noise, t_s);
    %%%%%%%% add controlller (oe_sys) & simlation
    n_ori.controllers = {};
    n_ori.add_controller(c_n, ss(d2c(sys_oe_set{itr},'foh')), Q, R);
    sys_ori_oe = n_ori.get_sys_controlled(sys_ori);
    y_oe_set(:, :, itr) = lsim(sys_ori_oe(ob_y_p, ID_in_p), sim_noise, t_s);
    xhat_oe_set(:, :, itr) = lsim(sys_ori_oe(ob_xhat_p, ID_in_p), sim_noise, t_s);
    %%%%%%%% add controlller (balred) & simlation
    n_ori.controllers = {};
    n_ori.add_controller(c_n, balred(sys_env, state), Q, R);
    sys_ori_ideal = n_ori.get_sys_controlled(sys_ori);
    y_ideal_set(:, :, itr) = lsim(sys_ori_ideal(ob_y_p, ID_in_p), sim_noise, t_s);
    xhat_ideal_set(:, :, itr) = lsim(sys_ori_ideal(ob_xhat_p, ID_in_p), sim_noise, t_s);
end
%% drawing

H1 = figure_config.set_figure_retro('theta');
yhat_ID_t = figure_config.plot_retro2(H1.axes, t_s, ...
                                        mean(y_ID_set(:,1,:),3), mean(xhat_ID_set(:,1,:),3), ...
                                        {'r-', 'linewidth', 1.0});
yhat_arx_t = figure_config.plot_retro2(H1.axes, t_s, ...
                                        mean(y_arx_set(:,1,:),3), mean(xhat_arx_set(:,1,:),3), ...
                                        {'b-', 'linewidth', 1.0});
yhat_armax_t = figure_config.plot_retro2(H1.axes, t_s, ...
                                        mean(y_armax_set(:,1,:),3), mean(xhat_armax_set(:,1,:),3), ...
                                        {'g-', 'linewidth', 1.0});
yhat_oe_t = figure_config.plot_retro2(H1.axes, t_s, ...
                                        mean(y_oe_set(:,1,:),3), mean(xhat_oe_set(:,1,:),3), ...
                                        {'k-', 'linewidth', 1.0});
yhat_ideal_t = figure_config.plot_retro2(H1.axes, t_s, ...
                                        mean(y_ideal_set(:,1,:),3), mean(xhat_ideal_set(:,1,:),3), ...
                                        {'m-', 'linewidth', 1.0});

H2 = figure_config.set_figure_retro('omega');
yhat_ID_o = figure_config.plot_retro2(H2.axes, t_s, ...
                                        mean(y_ID_set(:,2,:),3), mean(xhat_ID_set(:,2,:),3), ...
                                        {'r-', 'linewidth', 1.0});
yhat_arx_o = figure_config.plot_retro2(H2.axes, t_s, ...
                                        mean(y_arx_set(:,2,:),3), mean(xhat_arx_set(:,2,:),3), ...
                                        {'b-', 'linewidth', 1.0});
yhat_armax_o = figure_config.plot_retro2(H2.axes, t_s, ...
                                        mean(y_armax_set(:,2,:),3), mean(xhat_armax_set(:,2,:),3), ...
                                        {'g-', 'linewidth', 1.0});
yhat_oe_o = figure_config.plot_retro2(H2.axes, t_s, ...
                                        mean(y_oe_set(:,2,:),3), mean(xhat_oe_set(:,2,:),3), ...
                                        {'k-', 'linewidth', 1.0});
yhat_ideal_o = figure_config.plot_retro2(H2.axes, t_s, ...
                                        mean(y_ideal_set(:,2,:),3), mean(xhat_ideal_set(:,2,:),3), ...
                                        {'m-', 'linewidth', 1.0});                                    


H3 = figure_config.set_figure_bode('');
line_env = figure_config.plot_bode(H3.axes,omega,mag_env,pha_env);
line_ID = figure_config.plot_bode(H3.axes,omega,mean(Mag_ID,2),mean(Pha_ID,2),{'r-'});
line_arx = figure_config.plot_bode(H3.axes,omega,mean(Mag_arx,2),mean(Pha_arx,2),{'b-'});
line_armax = figure_config.plot_bode(H3.axes,omega,mean(Mag_armax,2),mean(Pha_armax,2),{'g-'});
line_oe = figure_config.plot_bode(H3.axes,omega,mean(Mag_oe,2),mean(Pha_oe,2),{'k-'});
[mag_ideal, pha_ideal] = bode_data(balred(sys_env,state),omega);
line_ideal = figure_config.plot_bode(H3.axes,omega,mag_ideal,pha_ideal,{'m-'});

legend(H3.axes.ax2,[line_env,line_ideal,line_arx,line_armax,line_oe,line_ID],...
        'original','balred','ARX','ARMAX','OE','ID')
% % % % % % % % save
% try
%     for itr = 1:4
%     figure(itr);
%     f = gcf;
%     savefig(f, strcat(root,'\',f.Name), 'compact')
%     end
%     close all
%     clear H1 H2 H3
%     save(strcat(root,'\data'))
% catch    
%     disp('Error !!')
% end
    
% end    
%% local 
function [sys1, sys2, sys3] = parfor_func(data,state)
    sys1 = armax(data, [state, state+1, state, 0]);
    sys2 = arx(data, [state, state+1, 0]);
    sys3 = oe(data, [state+1, state, 0]);
end

function [mag, pha, omega] = bode_data(sys,omega)
    if nargin < 2
        [mag, pha, omega] = bode(sys);
    else
        [mag, pha] = bode(sys, omega);
    end
    mag = squeeze(mag(1,1,:));
    pha = squeeze(pha(1,1,:));
end

function parfor_func_all()
   S = load('parfor_data.mat'); 
end
