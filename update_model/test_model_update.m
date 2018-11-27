clear
close all
%% Generate Network
seed = 10;
Node_number = 30;
n_ori = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
% n_ori.plot()
%% control & noise Node
c_n = 1;
n_n = [2:30];
%% signal power
id_in_p = 10;
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
add_io_node = unique([c_n,n_n]);
for nnn = add_io_node
    sys_ori = n_ori.add_io(sys_ori, nnn, strcat('node',num2str(nnn)));
end
[sys_local, sys_env] = n_ori.get_sys_local(c_n);
%% LQR params
Q = diag([1, 1]);
R = 1e-3;
%% Identificaiton model
state = 6;
in = 1;
out = 1;
%% Noise for compare response
Ts_s = 0.01;
Tf = 200;
t_s = 0:Ts_s:Tf;
rng(28);
sim_N = length(t_s);
sim_noise = randn(sim_N, 2);
%% ini
% res
N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;

rng('shuffle')

d = randn(N,2+numel(n_n)*2);
d(:,1:2) = d(:,1:2)*id_in_p*100;
d(:,3:end) = d(:,3:end)*noise_p;

n_ori.controllers = {};
n_ori.add_controller(c_n, Q, R);
sys_ori_con = n_ori.get_sys_controlled(sys_ori);

model_local = n_ori.controllers{1}.sys_K({'w'},{'v'});
Cz = zeros(1,order(model_local));
Cz(2) = 1;
model_spem_ini = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, in, out), model_local, Cz));
model_oe_ss_ini = model_ss(gen_ss_tridiag(state, in, out));
init_sys = model_spem_ini.gen_ss.gen_ss.get_sys();

model_spem_ini.max_iter = 5e3;
model_oe_ss_ini.max_iter = 5e3;

job1 = parfeval(@ID_func,1, sys_ori_con, d, t, ID_in_p, Noise, ob_v_p, ob_w_p, 'spem', model_spem_ini);
job2 = parfeval(@ID_func,1, sys_ori_con, d, t, ID_in_p, Noise, ob_v_p, ob_w_p, 'oe_ss', model_oe_ss_ini);
job3 = parfeval(@ID_func,2, sys_ori_con, d, t, ID_in_p, Noise, ob_v_p, ob_w_p, 'oe_tf', init_sys);
% job1 = parfeval(@model_spem.fit,1,t,[w_ori,v_ori],zeros(N, 1));
% job2 = parfeval(@model_oe_ss.fit,1,t,w_ori,v_ori);
% job3 = parfeval(@oe,1,iddata(v_ori,w_ori,Ts),init_sys);

while true
    if wait(job1) && wait(job2) && wait(job3)
        break;
    end
end

IDsys_ini_spem  = fetchOutputs(job1);
IDsys_ini_oe_ss = fetchOutputs(job2);
[IDsys_ini_oe_tf,IDsys_ini_oe_tf_id] = fetchOutputs(job3);

% ini spem
[y_ini_spem, ytilde_ini_spem, controller_ini_spem, sys_ini_spem] = sim_res(n_ori, sys_ori, c_n, IDsys_ini_spem, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, Q, R);
% ini oe_ss
[y_ini_oe_ss, ytilde_ini_oe_ss, controller_ini_oe_ss, sys_ini_oe_ss] = sim_res(n_ori, sys_ori, c_n, IDsys_ini_oe_ss, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, Q, R);
% ini oe_tf
[y_ini_oe_tf, ytilde_ini_oe_tf, controller_ini_oe_tf, sys_ini_oe_tf] = sim_res(n_ori, sys_ori, c_n, IDsys_ini_oe_tf, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, Q, R);

%% second
model_local = controller_ini_spem{1}.sys_K({'w'},{'v'});
Cz = zeros(1,order(model_local));
Cz(2) = 1;
model_spem_sec = model_ss(gen_ss_rectifier(gen_ss_tridiag(state, in, out), model_local, Cz));
model_oe_ss_sec = model_ss(gen_ss_tridiag(state, in, out));
model_spem_sec.gen_ss.gen_ss.set_sys(IDsys_ini_spem);
model_oe_ss_sec.gen_ss.set_sys(IDsys_ini_oe_ss);


model_spem_sec.max_iter = 5e3;
model_oe_ss_sec.max_iter = 5e3;

d = randn(N,2+numel(n_n)*2);
d(:,1:2) = d(:,1:2)*id_in_p;
d(:,3:end) = d(:,3:end)*noise_p;

job1 = parfeval(@ID_func,1, sys_ini_spem, d, t, ID_in_p, Noise, ob_v_p, ob_w_p, 'spem', model_spem_sec);
job2 = parfeval(@ID_func,1, sys_ini_oe_ss, d, t, ID_in_p, Noise, ob_v_p, ob_w_p, 'oe_ss', model_oe_ss_sec);
job3 = parfeval(@ID_func,2, sys_ini_oe_tf, d, t, ID_in_p, Noise, ob_v_p, ob_w_p, 'oe_tf', IDsys_ini_oe_tf_id);
% job1 = parfeval(@model_spem.fit,1,t,[w_ori,v_ori],zeros(N, 1));
% job2 = parfeval(@model_oe_ss.fit,1,t,w_ori,v_ori);
% job3 = parfeval(@oe,1,iddata(v_ori,w_ori,Ts),init_sys);

while true
    if wait(job1) && wait(job2) && wait(job3)
        break;
    end
end

IDsys_sec_spem  = fetchOutputs(job1);
IDsys_sec_oe_ss = fetchOutputs(job2);
[IDsys_sec_oe_tf, IDsys_sec_oe_tf_id] = fetchOutputs(job3);

% sec spem
[y_sec_spem, ytilde_sec_spem, controller_sec_spem, sys_sec_spem] = sim_res(n_ori, sys_ori, c_n, IDsys_sec_spem, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, controller_ini_spem{1}.K);
% sec oe_ss
[y_sec_oe_ss, ytilde_sec_oe_ss, controller_sec_oe_ss, sys_sec_oe_ss] = sim_res(n_ori, sys_ori, c_n, IDsys_sec_oe_ss, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, controller_ini_oe_ss{1}.K);
% sec oe_tf
[y_sec_oe_tf, ytilde_sec_oe_tf, controller_sec_oe_tf, sys_sec_oe_tf] = sim_res(n_ori, sys_ori, c_n, IDsys_sec_oe_tf, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, controller_ini_oe_tf{1}.K);
%% third
% sec spem
[y_thi_spem, ytilde_thi_spem, controller_thi_spem, sys_thi_spem] = sim_res(n_ori, sys_ori, c_n, IDsys_sec_spem, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, Q, R);
% sec oe_ss
[y_thi_oe_ss, ytilde_thi_oe_ss, controller_thi_oe_ss, sys_thi_oe_ss] = sim_res(n_ori, sys_ori, c_n, IDsys_sec_oe_ss, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, Q, R);
% sec oe_tf
[y_thi_oe_tf, ytilde_thi_oe_tf, controller_thi_oe_tf, sys_thi_oe_tf] = sim_res(n_ori, sys_ori, c_n, IDsys_sec_oe_tf, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, Q, R);

%% local 
function [y, ytilde, controller_info, sys_con] = sim_res(n_ori, sys_ori, c_n, sys_id, sim_noise, t_s, ob_y_p, ob_xhat_p, ID_in_p, Q, R)
    try
        n_ori.controllers = {};
        if nargin == 11
            n_ori.add_controller(c_n, sys_id, Q, R);
        else
            n_ori.add_controller2(c_n, sys_id, Q);
        end
        controller_info = n_ori.controllers;
        sys_con = n_ori.get_sys_controlled(sys_ori);
        y = lsim(sys_con(ob_y_p, ID_in_p), sim_noise, t_s);
        ytilde = lsim(sys_con(ob_xhat_p, ID_in_p), sim_noise, t_s);
    catch ME
       y = [];
       ytilde = [];
       controller_info = ME;
       sys_con = [];
    end
end

function [sys_out, sys_id] = ID_func(sys, d, t, ID_in_p, Noise, ob_v_p, ob_w_p, method, model)
    v = lsim(sys(ob_v_p, cat(2,ID_in_p,Noise)), d, t);
    w = lsim(sys(ob_w_p, cat(2,ID_in_p,Noise)), d, t);
    switch method
        case 'spem'
            model.fit(t, [w, v], zeros(size(w)));
            sys_out = model.gen_ss.gen_ss.get_sys();
        case 'oe_ss'
            model.fit(t, w, v);
            sys_out = model.gen_ss.get_sys(); 
        case 'oe_tf'
            Ts = t(2)-t(1);
            if model.Ts == 0
                model = c2d(model,Ts,'foh');
            end
            sys_id = oe(iddata(v,w,Ts), model); % model を初期システムとしている．
            sys_out = ss(d2c(sys_id));
    end
    if ~exist('sys_id','var')
        sys_id = [];
    end
end
