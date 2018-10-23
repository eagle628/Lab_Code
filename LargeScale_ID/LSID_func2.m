%% NOTE
% ARX, ARMAX of System Idnetification Toolbox cannnot use in PARFOR LOOP.
% So, the part do local function
% SS-OE 比較
%% 
function LSID_func2(Node_number, noise_p, n_n, root, max_iter, controller_number, c_n)
% clear 
close all
% Node_number=3; noise_p = 0.1; n_n=[2,3]; root=''; max_iter=4; controller_number=[3]; c_n=1;
%% Generate Network
% seed = 3;% node numver 3
seed = 10;% node numver 30
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
for itr = 1:numel(controller_number) % node1 is target local system
    n_ori.add_controller( controller_number(itr), Q, R);
end
sys_ori_c1 = n_ori.get_sys_controlled(sys_ori);
%% Generate v & w
N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;
%% fitting & simlation
%%%%%%%%%% fitting config
% max_iter = 20;
state = 4;
in = 1;
out = 1;

%%%% Extend response simlation
Ts_s = 0.01;
Tf = 200;
t_s = 0:Ts_s:Tf;
rng(28);
sim_N = length(t_s);
sim_noise = randn(sim_N, 2);
%%%%%%%% controller parameter
Q = kron(eye(1*numel(1)),diag([1,1000]));
R = kron(eye(1*numel(1)),diag([1e-3]));


parfor_progress(max_iter);
%%%%%% init sys (coloum vector)
rng('shuffle')

save('parfor_data');

% get memory
Mag_ID = zeros(length(omega), max_iter);
Pha_ID = zeros(length(omega), max_iter);
Mag_arx = zeros(length(omega), max_iter);
Pha_arx = zeros(length(omega), max_iter);
Mag_armax = zeros(length(omega), max_iter);
Pha_armax = zeros(length(omega), max_iter);
Mag_oe = zeros(length(omega), max_iter);
Pha_oe = zeros(length(omega), max_iter);
Mag_oe_ss = zeros(length(omega), max_iter);
Pha_oe_ss = zeros(length(omega), max_iter);

sys_arx_set = cell(max_iter, 1);
sys_armax_set = cell(max_iter, 1);
sys_oe_set = cell(max_iter, 1);
sys_ID_set = cell(max_iter, 1);
sys_oe_ss_set = cell(max_iter, 1);

best_cost_set = zeros(max_iter, 1);
ID_cost_set = zeros(max_iter, 1);
oe_ss_cost_set = zeros(max_iter, 1);

% get memory
y_ID_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
y_arx_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
y_armax_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
y_oe_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
y_ideal_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
y_ss_oe_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_ID_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_arx_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_armax_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_oe_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_ideal_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);
xhat_ss_oe_set = zeros(sim_N, 2*numel(ob_y_p), max_iter);

parfor itr = 1:max_iter
    [ best_cost_set(itr), sys_armax_set{itr}, sys_arx_set{itr}, sys_oe_set{itr},...
            Mag_armax(:, itr), Pha_armax(:, itr), Mag_arx(:, itr), Pha_arx(:, itr), Mag_oe(:, itr), Pha_oe(:, itr),...
            sys_ID_set{itr}, ID_cost_set(itr), Mag_ID(:, itr), Pha_ID(:, itr),...
            sys_oe_ss_set{itr}, oe_ss_cost_set(itr), Mag_oe_ss(:, itr), Pha_oe_ss(:, itr),...
            y_ID_set(:,:,itr), xhat_ID_set(:,:,itr), y_armax_set(:,:,itr), xhat_armax_set(:,:,itr), y_arx_set(:,:,itr),...
            xhat_arx_set(:,:,itr), y_oe_set(:,:,itr), xhat_oe_set(:,:,itr), y_ideal_set(:,:,itr), xhat_ideal_set(:,:,itr), y_ss_oe_set(:,:,itr), xhat_ss_oe_set(:,:,itr)] = parfor_func_all()
    
    parfor_progress();
end
parfor_progress(0);

%% simlation
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
yhat_ss_oe_t = figure_config.plot_retro2(H1.axes, t_s, ...
                                        mean(y_ss_oe_set(:,1,:),3), mean(xhat_ss_oe_set(:,1,:),3), ...
                                        {'c-', 'linewidth', 1.0});

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
yhat_ss_oe_o = figure_config.plot_retro2(H2.axes, t_s, ...
                                        mean(y_ss_oe_set(:,2,:),3), mean(xhat_ss_oe_set(:,2,:),3), ...
                                        {'c-', 'linewidth', 1.0});                                    


H3 = figure_config.set_figure_bode('');
line_env = figure_config.plot_bode(H3.axes,omega,mag_env,pha_env);
line_ID = figure_config.plot_bode(H3.axes,omega,mean(Mag_ID,2),mean(Pha_ID,2),{'r-'});
line_arx = figure_config.plot_bode(H3.axes,omega,mean(Mag_arx,2),mean(Pha_arx,2),{'b-'});
line_armax = figure_config.plot_bode(H3.axes,omega,mean(Mag_armax,2),mean(Pha_armax,2),{'g-'});
line_oe = figure_config.plot_bode(H3.axes,omega,mean(Mag_oe,2),mean(Pha_oe,2),{'k-'});
[mag_ideal, pha_ideal] = bode_data(balred(sys_env,state),omega);
line_ideal = figure_config.plot_bode(H3.axes,omega,mag_ideal,pha_ideal,{'m-'});
line_oe_ss = figure_config.plot_bode(H3.axes,omega,mean(Mag_oe_ss,2),mean(Pha_oe_ss,2),{'c-'});

legend(H3.axes.ax2,[line_env,line_ideal,line_arx,line_armax,line_oe,line_ID,line_oe_ss],...
        'original','balred','ARX','ARMAX','OE-TF','ID','OE-SS')
%% save

try
    for itr = 1:4
        figure(itr);
        f = gcf;
        savefig(f, strcat(root,'\',f.Name), 'compact')
    end
    close all
    clear H1 H2 H3
    save(strcat(root,'\data'))
catch    
    disp('Error !!')
end
    
end    


%% local 

function [mag, pha, omega] = bode_data(sys,omega)
    if nargin < 2
        [mag, pha, omega] = bode(sys);
    else
        [mag, pha] = bode(sys, omega);
    end
    mag = squeeze(mag(1,1,:));
    pha = squeeze(pha(1,1,:));
end

function [ best_cost, sys_tf_armax, sys_tf_arx, sys_tf_oe,...
            mag_armax, pha_armax, mag_arx, pha_arx, mag_oe, pha_oe,...
            sys_ID, ID_cost, mag_ID, pha_ID,...
            sys_ss_oe, oe_ss_cost, mag_ss_oe, pha_ss_oe,...
            y_ID, xhat_ID, y_armax, xhat_armax, y_arx,...
            xhat_arx, y_oe, xhat_oe, y_ideal, xhat_ideal, y_ss_oe, xhat_ss_oe  ] = parfor_func_all()
        
    % load data
    S = load('parfor_data.mat');
    
    % start fitting
   % Response of v&w 
    rng('shuffle');
    d = randn(S.N,2+numel(S.n_n)*2);
    d(:,1:2) = (d(:,1:2)+0.5)*S.id_in_p;
    d(:,3:end) = d(:,3:end)*S.noise_p;
  
    v = lsim(S.sys_ori_c1(S.ob_v_p, cat(2,S.ID_in_p,S.Noise)), d, S.t);
    w = lsim(S.sys_ori_c1(S.ob_w_p, cat(2,S.ID_in_p,S.Noise)), d, S.t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % best cost calculate
    m_Best = model_ss(gen_ss_rectifier(gen_ss_all(S.state, S.in, S.out), S.sys_local_vw, [0 1]));
    m_Best.set_sys(balred(S.sys_env ,S.state));
    best_cost = m_Best.eval_func(S.t, [w, v], zeros(S.N, 1));
%     fprintf('Best Cost is %e.', best_cost);
%     fprintf('\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% Choice ID SS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For Identification
%     m = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out), sys_local_vw, [0 1]));
%     m = model_ss(gen_ss_rectifier(gen_ss_canonical(state, in, out), sys_local_vw, [0 1]));
    m = model_ss(gen_ss_rectifier(gen_ss_tridiag(S.state, S.in, S.out), S.sys_local_vw, [0 1]));
%     m = model_ss(gen_ss_rectifier(gen_ss_canonical_zero(state, in, out), sys_local_vw, [0 1]));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%w%%%%%%%%%%%%%%
    % Open Loop ID
    data = iddata(v, w, S.Ts);
    init_sys = m.gen_ss.gen_ss.get_sys();
    sys_tf_armax = armax(data, c2d(init_sys, S.Ts, 'foh'));
    sys_tf_arx = arx(data, [S.state, S.state+1, 0]);
    sys_tf_oe = oe(data, c2d(init_sys, S.Ts, 'foh'));
    [mag_armax, pha_armax] = bode_data(sys_tf_armax,S.omega);
    [mag_arx, pha_arx] = bode_data(sys_tf_arx,S.omega);
    [mag_oe, pha_oe] = bode_data(sys_tf_oe,S.omega);
%     
%     m.set_sys(balred(sys_env,state));
%     m.set_sys(balred(ss(d2c(sys_arx ,'foh')),state));
%     theta = m.get_params_all();
    
%%%%%%%%%%%%%%% fitting
    iter_toler = 1e4;

    m.max_iter = iter_toler;
    m.str_display = 'none';
    % optimizer
    [theta, J] = m.fit(S.t, [w, v], zeros(S.N, 1));
    sys_ID = m.gen_ss.gen_ss.get_sys();
    ID_cost = m.eval_func(S.t, [w, v], zeros(S.N, 1));
    % model character
%     [mag, pha] = bode(model, omega);
%     mag = squeeze(mag(1,1,:));
%     pha = squeeze(pha(1,1,:));
%     Mag_ID(:, itr) = mag;
%     Pha_ID(:, itr) = pha;
    [mag_ID, pha_ID] = bode_data(sys_ID, S.omega);
    
    m1 = model_ss(gen_ss_tridiag(S.state, S.in, S.out));
    m1.max_iter = iter_toler;
    m1.str_display = 'none';
    m1.fit(S.t, w, v);
    sys_ss_oe = m1.gen_ss.get_sys();
    [mag_ss_oe, pha_ss_oe] = bode_data(sys_ss_oe, S.omega);
    oe_ss_cost = m1.eval_func(S.t, w, v);
    
    % end fitting 
    
    % start compare response
    %%%%%%%% add controller (ID) & simlation
    S.n_ori.controllers = {};
    S.n_ori.add_controller( S.c_n, sys_ID, S.Q, S.R);
    sys_ori_ID = S.n_ori.get_sys_controlled(S.sys_ori);
    y_ID = lsim(sys_ori_ID(S.ob_y_p, S.ID_in_p), S.sim_noise, S.t_s);
    xhat_ID = lsim(sys_ori_ID(S.ob_xhat_p, S.ID_in_p), S.sim_noise, S.t_s);
    %%%%%%%% add controlller (armax_sys) & simlation
    S.n_ori.controllers = {};
    S.n_ori.add_controller(S.c_n, ss(d2c(sys_tf_armax,'foh')), S.Q, S.R);
    sys_ori_armax = S.n_ori.get_sys_controlled(S.sys_ori);
    y_armax = lsim(sys_ori_armax(S.ob_y_p, S.ID_in_p), S.sim_noise, S.t_s);
    xhat_armax = lsim(sys_ori_armax(S.ob_xhat_p, S.ID_in_p), S.sim_noise, S.t_s);
    %%%%%%%% add controlller (arx_sys) & simlation
    S.n_ori.controllers = {};
    S.n_ori.add_controller(S.c_n, ss(d2c(sys_tf_arx,'foh')), S.Q, S.R);
    sys_ori_arx = S.n_ori.get_sys_controlled(S.sys_ori);
    y_arx = lsim(sys_ori_arx(S.ob_y_p, S.ID_in_p), S.sim_noise, S.t_s);
    xhat_arx = lsim(sys_ori_arx(S.ob_xhat_p, S.ID_in_p), S.sim_noise, S.t_s);
    %%%%%%%% add controlller (oe_sys) & simlation
    S.n_ori.controllers = {};
    S.n_ori.add_controller(S.c_n, ss(d2c(sys_tf_oe,'foh')), S.Q, S.R);
    sys_ori_oe = S.n_ori.get_sys_controlled(S.sys_ori);
    y_oe = lsim(sys_ori_oe(S.ob_y_p, S.ID_in_p), S.sim_noise, S.t_s);
    xhat_oe = lsim(sys_ori_oe(S.ob_xhat_p, S.ID_in_p), S.sim_noise, S.t_s);
    %%%%%%%% add controlller (balred) & simlation
    S.n_ori.controllers = {};
    S.n_ori.add_controller(S.c_n, balred(S.sys_env, S.state), S.Q, S.R);
    sys_ori_ideal = S.n_ori.get_sys_controlled(S.sys_ori);
    y_ideal = lsim(sys_ori_ideal(S.ob_y_p, S.ID_in_p), S.sim_noise, S.t_s);
    xhat_ideal = lsim(sys_ori_ideal(S.ob_xhat_p, S.ID_in_p), S.sim_noise, S.t_s);
    %%%%%%%% add controlller (ss_oe) & simlation
    S.n_ori.controllers = {};
    S.n_ori.add_controller(S.c_n, sys_ss_oe, S.Q, S.R);
    sys_ori_ss_oe = S.n_ori.get_sys_controlled(S.sys_ori);
    y_ss_oe = lsim(sys_ori_ss_oe(S.ob_y_p, S.ID_in_p), S.sim_noise, S.t_s);
    xhat_ss_oe = lsim(sys_ori_ss_oe(S.ob_xhat_p, S.ID_in_p), S.sim_noise, S.t_s);
    
end
