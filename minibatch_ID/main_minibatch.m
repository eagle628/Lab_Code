clear 
close all
%% Generate Network
seed = 3;
Node_number = 3;
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
n_ori.add_controller( c_n, Q, R);
sys_ori_c1 = n_ori.get_sys_controlled(sys_ori);
%% Generate v & w
N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;
%% minibatch
maxitr = 1;
state = 2*(Node_number-1);
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
    % best cost calculate
    Best_cost = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out), sys_local_vw));
    best_cost = Best_cost.eval_func(t, [w, v], zeros(N, 1), sys2params_all(sys_env, state, in, out));
    fprintf('Best Cost is %e.', best_cost);
    fprintf('\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % For Identification
    m = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out), sys_local_vw));
%     m = model_ss(gen_ss_rectifier(gen_ss_canonical(state, in, out), sys_local_vw));
%     m = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, in, out), sys_local_vw));
%     m = model_ss(gen_ss_rectifier(gen_ss_canonical_zero(state, in, out), sys_local_vw));
    % init
%     init_sys = sys_env;
    init_sys = d2c(arx(iddata(v,w,Ts), [state, state+1, 0])); 
%     init_sys = d2c(armax(iddata(v,w,Ts), [state, state+1, state, 0])); 
%     init_sys = d2c(oe(iddata(v,w,Ts), [state+1, state, 0])); 
%     init_params = sys2params_tri(init_sys, state, in, out); 
%     init_params = ones(3*state-2+state*(in+out)+in*out, 1)*1e-3;
%     m.add_fixed_params('theta_A_6', 0); % model.m
%     m.add_fixed_params('theta_A_9', 0); % model.m
%     init_params = nonzeros(init_params);
%     init_params = sys2params_canonical(ss(init_sys), state, in, out);
%     init_sys = rss(state);
%     init_params = sys2params_all(ss(init_sys), state, in, out);
%     init_params = sys2params_all(sys_env, state, in, out);
    init_params = sys2params_all(init_sys, state, in, out);
%     init_params = sys2params_canonical(init_sys, state, in, out);
%     init_params = sys2params_canonical_zero(init_sys, state, in, out);
    
%     m.add_fixed_params('theta_C_4', 0);
%%%%%%%%%%%%%%% fitting
    m.max_iter = 10000;
    % optimizer
    rng('shuffle')
    load('params.mat');
%     load('tekitou.mat')
%     [num,den] = tfdata(sys_env,'v');
%     den(2:end) = den(2:end)+randn(1,4)*0.1;
%     num = num+randn(1,5)*0.1;
%     sys = tf(num,den);
%     figure('position',[-600, 0, 600, 600])
%     bode(sys)
%     init_params = sys2params_all(ss(init_sys), state, in, out);
%     load('sys.mat')
%     while true
%         T = randn(state);
%         if (abs(det(T)) > 1e-6) && (abs(det(T)) < 1e6)
%             break;
%         end
%     end
%     init_params = sys2params_all(ss2ss(target_system, T), state, in, out);
%     init_params = sys2params_all(balreal(target_system), state, in, out);
%     theta= init_params;

%     init_params = ones(3*state-2+state*(in+out)+in*out, 1)*randn(1)*1e-4;init_params(end) = 1e-3;
%     for itr1 = 1: 10
%         [theta ,J] = m.fit_adam(t, [w, v], zeros(N, 1), theta, 1e-5, 0.9, 0.999, 1e-4, 0.8);
%         [theta ,J] = m.fit_adamax(t, [w, v], zeros(N, 1), theta, 1e-12, 0.9, 0.999, 0.8);
%         [theta ,J] = m.fit_rmsprop(t, [w, v], zeros(N, 1), theta, [], [], [], 0.8);
%         [theta ,J] = m.fit_SMORMS3(t, [w, v], zeros(N,1), theta, 1e-10, [], 0.8);
%         [theta, J] = m.fit_LMA(t, [w, v], zeros(N, 1), theta, [], []);
%         model1 = m.gen_ss.gen_ss.get_sys();
%         bode(sys_env,init_sys,model1)
%     end
%     [re_theta ,J] = m.fit_adamax(t, [w, v], zeros(N, 1), init_params, [], [], [], 0.8);
%     theta = 1e-4;
%     theta = init_params;
%     for itr1 = 1 : 2
        [theta,J] = m.fit_Santa_S(t, [w, v], zeros(N, 1), theta, 1e-10, [], [], [], @(t)t^2, 0.8);
%         model1 = m.gen_ss.gen_ss.get_sys();
%         bode(sys_env,init_sys,model1)
%     end
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

%% simlation
Ts_s = 0.01;
Tf = 200;
t_s = 0:Ts_s:Tf;
sim_noise = randn(length(t_s), 2);
% id
sys_ARMAX = d2c(armax(iddata(v,w,Ts), [state, state+1, state, 0])); 
%% controller parameter
Q = kron(eye(1*numel(1)),diag([1,1000]));
R = kron(eye(1*numel(1)),diag([1e-3]));
%% add controller (ID) & simlation
n_ori.controllers = {};
n_ori.add_controller( c_n, model1, Q, R);
sys_ori_ID = n_ori.get_sys_controlled(sys_ori);
y_ID = lsim(sys_ori_ID(ob_y_p, ID_in_p), sim_noise, t_s);
xhat_ID = lsim(sys_ori_ID(ob_xhat_p, ID_in_p), sim_noise, t_s);
%% add controlller (init_sys) & simlation
n_ori.controllers = {};
n_ori.add_controller(c_n, ss(init_sys), Q, R);
sys_ori_init = n_ori.get_sys_controlled(sys_ori);
y_init = lsim(sys_ori_init(ob_y_p, ID_in_p), sim_noise, t_s);
xhat_init = lsim(sys_ori_ID(ob_xhat_p, ID_in_p), sim_noise, t_s);
%% add controlller (init_sys) & simlation
n_ori.controllers = {};
n_ori.add_controller(c_n, ss(sys_ARMAX), Q, R);
sys_ori_armax = n_ori.get_sys_controlled(sys_ori);
y_armax = lsim(sys_ori_armax(ob_y_p, ID_in_p), sim_noise, t_s);
xhat_armax = lsim(sys_ori_armax(ob_xhat_p, ID_in_p), sim_noise, t_s);
%% drawing
H1 = figure_config.set_figure_retro('theta');
ytilde_ID = figure_config.plot_retro2(H1.axes, t_s, y_ID(:,1), xhat_ID(:,1), {'r-', 'linewidth', 2.0});
ytilde_init = figure_config.plot_retro2(H1.axes, t_s, y_init(:,1), xhat_init(:,1), {'b:', 'linewidth', 2.0});
ytilde_armax = figure_config.plot_retro2(H1.axes, t_s, y_armax(:,1), xhat_armax(:,1), {'g--', 'linewidth', 1.6});

H2 = figure_config.set_figure_retro('omega');
ytilde_ID = figure_config.plot_retro2(H2.axes, t_s, y_ID(:,2), xhat_ID(:,2), {'r-', 'linewidth', 2.0});
ytilde_init = figure_config.plot_retro2(H2.axes, t_s, y_init(:,2), xhat_init(:,2), {'b:', 'linewidth', 2.0});
ytilde_armax = figure_config.plot_retro2(H2.axes, t_s, y_armax(:,2), xhat_armax(:,2), {'g--', 'linewidth', 1.6});

figure('Name','Bode')
bode(sys_env,init_sys,model1,sys_ARMAX)
legend('original','Init','ID','ARMAX')

figure('Name','Cost History')
plot(nonzeros(J))
%% local function
function init_params = sys2params_tri(init_sys, state, in, out)
    init_sys = canon(balred(ss(init_sys),state),'modal');
    init_params = zeros(3*state-2+state*(in+out)+in*out, 1);
    init_params(1:3*state-2) = [diag(init_sys.A)', diag(init_sys.A, -1)', diag(init_sys.A, 1)'];
    init_params(3*state-2+1:3*state-2+state*in) = reshape(init_sys.B, state*in, 1);
    init_params(3*state-2+state*in+1:3*state-2+state*in+out*state) = reshape(init_sys.C, out*state, 1);
    init_params(3*state-2+state*in+out*state+1:3*state-2+state*in+out*state+in*out) = reshape(init_sys.D, out*in, 1);
end

function init_params = sys2params_all(init_sys, state, in, out)
    init_sys = balred(ss(init_sys),state);
    init_params = zeros(state^2+state*(in+out)+in*out, 1);
    init_params(1:state^2) = reshape(init_sys.A, state^2, 1);
    init_params(state^2+1:state^2+state*in) = reshape(init_sys.B, state*in, 1);
    init_params(state^2+state*in+1:state^2+state*in+out*state) = reshape(init_sys.C, out*state, 1);
    init_params(state^2+state*in+out*state+1:state^2+state*in+out*state+out*in) = reshape(init_sys.D, out*in, 1);    
end

function init_params = sys2params_canonical(init_sys, state, in, out)
    init_sys = balred(ss(init_sys),state);
%     init_sys = canon(init_sys, 'companion');
    init_sys = ctrbcanon(init_sys);
    init_params = zeros(state + state*((in-1)+out) + in*out, 1);
    init_params(1:state) = init_sys.A(end, :);
    init_params(state+1:state+state*(in-1)) = reshape(init_sys.B(:,2:end), state*(in-1), 1);
    init_params(state+state*(in-1)+1:state+state*(in-1)+out*state) = reshape(init_sys.C, out*state, 1);
    init_params(state+state*(in-1)+out*state+1:state+state*(in-1)+out*state+out*in) = reshape(init_sys.D, out*in, 1);    
end

function init_params = sys2params_canonical_zero(init_sys, state, in, out)
    init_sys = balred(ss(init_sys),state);
%     init_sys = canon(init_sys, 'companion');
    init_sys = ctrbcanon(init_sys);
    init_params = zeros(state + state*((in-1)+out) + in*out, 1);
    init_params(1:state) = init_sys.A(end, :);
    init_params(state+1:state+state*(in-1)) = reshape(init_sys.B(:,2:end), state*(in-1), 1);
    init_params(state+state*(in-1)+1:state+state*(in-1)+out*state) = reshape(init_sys.C, out*state, 1);
    init_params(state+state*(in-1)+out*state+1:state+state*(in-1)+out*state+out*in) = reshape(init_sys.D, out*in, 1);    
    init_params(state+1) = [];
end

function csys = ctrbcanon(sys)
    % Only SISO(for continuous)
    [num,den] = tfdata(sys ,'v');
    dim = order(sys);
    A = [zeros(dim-1,1), eye(dim-1); -fliplr(den(2:end))];
    B = [zeros(dim-1,1); 1];
    C = fliplr(num(2:end) - den(2:end)*num(1));
    D = num(1);
    csys = ss(A,B,C,D);
end
%%
% m.add_fixed_params('theta_D_1', 0);
