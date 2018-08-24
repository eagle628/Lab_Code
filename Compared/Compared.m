% Environment Identification
%{
Compare Response
Origianl System Response
Simple Retrofit Control
Extend Retrofit Control
%}
%% initiallize workspace
clear 
close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 30;
seed = 10;
%rng('shuffle');
%seed = randi(1000,1,1);
m = [1,2];
d = [2,10]*1e-2;
b = 1;
Y = [ 1, 5];
% Y = [ 0, 1];
r_c = 0.1;
n_ori = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% Controlled node number
c_n = 1;
%% simulation noise point & power
noise_point = c_n;
noise_power = 1;
%% Node Name
sim_noise_p = cell(1,numel(noise_point));
for i = 1 :  numel(noise_point)
    sim_noise_p(i) ={strcat('d_node',num2str(noise_point(i)))};
end
ob_y_p  = {strcat('y_node',num2str(c_n))};
ob_xhat_p = {'xhat_controlled1'};
ob_v_p = {strcat('v_node',num2str(c_n))};
%% add I/O port for original system
sys_ori = n_ori.get_sys();
add_io = unique([ c_n, noise_point]);
for i =  1 : numel(add_io)
    sys_ori = n_ori.add_io(sys_ori,add_io(i), strcat('node',num2str(add_io(i))));
end
%% Environment Character (Original)
% Chose Range Fruquency
min = -6;
max =  6;
[~, sys_env_h] = n_ori.get_sys_local(c_n);
sys_env = balred(sys_env_h,6);
[mag_env,pha_env,omega] = bode(sys_env_h,{10^min,10^max});
mag_env=squeeze(mag_env(1,1,:));
pha_env=squeeze(pha_env(1,1,:));
%% LQR Controller Parameter
Q = kron(eye(1*numel(c_n)),diag([1,1000]));
R = kron(eye(1*numel(c_n)),diag([1e-3]));
%% Add Simple RetroFit Controlelr
n_ori.controllers = {};
n_ori.add_controller( c_n, Q, R);
sys_simple = n_ori.get_sys_controlled(sys_ori);
sys_simple_K = n_ori.controllers{1}.sys_K;
%% Add Extend RetroFit Controlelr
n_ori.controllers = {};
n_ori.add_controller( c_n, sys_env_h, Q, R);
sys_extend = n_ori.get_sys_controlled(sys_ori);
sys_extend_K = n_ori.controllers{1}.sys_K;
%% Simlulation point
simulation_time = 500;%simulation_time
Ts_s = 0.01;
t_s = 0:Ts_s:simulation_time-Ts_s;
t_s = t_s';

sim_seed = 1024;
cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',length(t_s),'NumChannels',2*numel(noise_point),'RandomStream','mt19937ar with seed','Seed',sim_seed);
sim_noise = cn3()*noise_power;
clear cn3 

pattern = 1;
switch pattern
    case 1
        % original
        y_ori = lsim( sys_ori( ob_y_p, sim_noise_p), sim_noise, t_s);
        v_ori = lsim( sys_ori( ob_v_p, sim_noise_p), sim_noise, t_s);
        % Simple        
        y_simple = lsim( sys_simple( ob_y_p, sim_noise_p), sim_noise, t_s);
        yhat_simple = lsim( sys_simple( ob_xhat_p, sim_noise_p), sim_noise, t_s);
        v_simple = lsim( sys_simple( ob_v_p, sim_noise_p), sim_noise, t_s);
        % Extend
        y_extend = lsim( sys_extend( ob_y_p, sim_noise_p), sim_noise, t_s);
        yhat_extend = lsim( sys_extend( ob_xhat_p, sim_noise_p), sim_noise, t_s);
        v_extend = lsim( sys_extend( ob_v_p, sim_noise_p), sim_noise, t_s);
    case 2
        % original
        y_ori = impulse( sys_ori( ob_y_p, sim_noise_p), t_s);
        yhat_ori = impulse( sys_ori( 'x', sim_noise_p), t_s);
        v_ori = impulse( sys_ori( ob_v_p, sim_noise_p), t_s);
        % Simple
        y_simple = impulse( sys_simple( ob_y_p, sim_noise_p), t_s);
        yhat_simple = impulse( sys_simple( ob_xhat_p, sim_noise_p), t_s);
        v_simple = impulse( sys_simple( ob_v_p, sim_noise_p), t_s);
        %extend
        y_extend = impulse( sys_extend( ob_y_p, sim_noise_p), t_s);
        yhat_extend = impulse( sys_extend( ob_xhat_p, sim_noise_p), t_s);
        v_extend = impulse( sys_extend( ob_v_p, sim_noise_p), t_s);
    otherwise
        disp('Error')
end

%% Drawing
state = 2;
check = 1;
fig3 = figure_config.set_figure_retro('Response');

% Original
figure_config.plot_retro2( fig3.axes, t_s, y_ori(:,state,check), [], {'g-','linewidth',3.0});
% Simple
figure_config.plot_retro2( fig3.axes, t_s, y_simple(:,state,check), yhat_simple(:,state,check), {'r-','linewidth',2.5});
% Extend
figure_config.plot_retro2( fig3.axes, t_s, y_extend(:,state,check), yhat_extend(:,state,check), {'b-','linewidth',2.0});

ori = fig3.axes.ax1.Children(end);
simple = fig3.axes.ax1.Children(end-1);
extend = fig3.axes.ax1.Children(end-2);
legend(fig3.axes.ax1,[ori,simple,extend],'original','simple','extend')