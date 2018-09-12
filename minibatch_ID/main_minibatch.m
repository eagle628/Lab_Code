clear 
close all
%% Generate Network
seed = 3;
Node_number = 2;
n_ori = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% control & noise Node
c_n = 1;
n_n = [2];
%% signal power
id_in_p = 1;
noise_p = 0;
%% Node Name
ob_v  = {strcat('v_node',num2str(c_n))};
ob_w  = {strcat('w_node',num2str(c_n))};
ID_in = {strcat('d_node',num2str(c_n))};
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
% Q = kron(eye(1*numel(1)),diag([1,1000]));
% R = kron(eye(1*numel(1)),diag([1e-3]));
% n_ori.add_controller( 2, Q, R);
% sys_ori = n_ori.get_sys_controlled(sys_ori);
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

%
for itr = 1:maxitr
    d = randn(N,2+numel(n_n)*2);
    d(:,1:2) = (d(:,1:2)+0.5)*id_in_p;
    d(:,3:end) = d(:,3:end)*noise_p;
    % Response of v&w 
    rng(28);
    v = lsim(sys_ori(ob_v, cat(2,ID_in,Noise)), d, t);
    w = lsim(sys_ori(ob_w, cat(2,ID_in,Noise)), d, t);

    m = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, in, out), sys_local_vw));
%     m = model_ss(gen_ss_rectifier(gen_ss_all(state, in, out), sys_local_vw));
%     m = model_ss(gen_ss_rectifier(gen_ss_canonical(state, in, out), sys_local_vw));
    m.max_iter = 5e3;
    %%%%%%%%%%%%%%%%%%%%%
    best_cost = m.eval_func(t, [w, v], zeros(N, 1), sys2params_tri(sys_env, state, in, out));
%     best_cost = m.eval_func(t, [w, v], zeros(N, 1), sys2params_all(sys_env, state, in, out));
%     best_cost = m.eval_func(t, [w, v], zeros(N, 1), sys2params_canonical(sys_env, state, in, out));
    fprintf('Best Cost is %e.', best_cost);
    fprintf('\n');
    %%%%%%%%%%%%%%%%%%%%%
    % init
    init_sys = d2c(arx(iddata(v,w,Ts), [state, state+1, 0])); 
    init_params = sys2params_tri(init_sys, state, in, out); 
%     init_params = ones(3*state-2+state*(in+out)+in*out, 1)*1e-3;
%     m.add_fixed_params('theta_A_6', 0); % model.m
%     m.add_fixed_params('theta_A_9', 0); % model.m
%     init_params = nonzeros(init_params);
%     init_params = sys2params_canonical(ss(d2c((init_sys))), state, in, out);
%     init_sys = rss(state);
%     init_params = sys2params_all(ss(d2c(init_sys)), state, in, out);
%     init_params = sys2params_all(init_sys, state, in, out);
    
    % optimizer
    [theta ,J] = m.fit_adam(t, [w, v], zeros(N, 1), init_params, 1e-5, 0.9, 0.999, [], 0.5, 1);
%     [theta ,J] = m.fit_adamax(t, [w, v], zeros(N, 1), init_params, [], [], [], 0.6);
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
    
    [mag, pha] = bode(init_sys, omega);
    mag = squeeze(mag(1,1,:));
    pha = squeeze(pha(1,1,:));
figure_config.plot_bode(H.axes,omega,mag,pha,{'g:','linewidth',2.0});

figure_config.plot_bode(H.axes,omega,Mag,Pha,{'b-','linewidth',0.8});
    
figure('Name','Cost History')
semilogy(J)

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
    init_sys = canon(init_sys, 'companion');
    init_params = zeros(state + state*((in-1)+out) + in*out, 1);
    init_params(1:state) = init_sys.A(:, end);
    init_params(state+1:state+state*(in-1)) = reshape(init_sys.B(:,2:end), state*(in-1), 1);
    init_params(state+state*(in-1)+1:state+state*(in-1)+out*state) = reshape(init_sys.C, out*state, 1);
    init_params(state+state*(in-1)+out*state+1:state+state*(in-1)+out*state+out*in) = reshape(init_sys.D, out*in, 1);    
end

%%
%     m.add_fixed_params('theta_D_1', 0); % model.m
