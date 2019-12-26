close all
clear
%% condition
data.node_number = 2;
data.net_seed = 6;
data.net_local = 1;
data.Ts = 0.1;

data.controller_n = 4;

data.eval_seed = 6;

data.ga_seed = 12;
data.pso_seed = 12;
data.train_input_seed = 28;
data.controller_struct_enable = true;
data.train_input_form = true;
% data.controller_gain_enable = true;

SSS = load('2019-12-09=00-05-20');
data = SSS.data;


%% define model
% % % network version
net = network_swing_simple(data.node_number, [1,2], [2,12]*1e-3, 1, [0.1,6], 0.8, data.net_seed);
model = swing_network_model(net, data.net_local, data.Ts);
% SLQR
net.add_controller(data.net_local, diag([1,1]), 1);
sys_all_slqr = net.get_sys();
sys_all_slqr = net.add_io(sys_all_slqr, data.net_local, 'node1');
sys_all_slqr = net.get_sys_controlled(sys_all_slqr);
slqr_K = net.controllers{1}.K;
net.controllers = {};
% ELQR
[~, sys_env] = net.get_sys_local(data.net_local);
env_n = data.controller_n - 2;
net.add_controller(data.net_local, balred(sys_env, env_n), diag([1,1]), 1);
sys_all_elqr = net.get_sys();
sys_all_elqr = net.add_io(sys_all_elqr, data.net_local, 'node1');
sys_all_elqr = net.get_sys_controlled(sys_all_elqr);
net.controllers = {};
%% define controller model
% set belief_N
belief_N = 1;
env_n = data.controller_n - 2;

if data.controller_struct_enable
    env_model = gen_ss_tridiag(env_n, 1, 1);
%     if data.controller_gain_enable
%         ss_model = gen_ss_Eretro_LQR_controller(env_model, c2d(model.sys_local, model.Ts, model.discrete_type));
%     else
    ss_model = gen_ss_Eretro_controller(env_model, c2d(model.sys_local, model.Ts, model.discrete_type));
%     end
else
    ss_model = gen_ss_tridiag(data.controller_n, model.ny*belief_N, model.nu);
end
%% Constructor
interpreter_model = observation_accumulater(model.ny, belief_N);
Target = Global_Target(model, ss_model, interpreter_model);

%%
multi_startpt_number = 1;
theta_set = zeros(multi_startpt_number, ss_model.N);

rng(data.pso_seed)
nvar = ss_model.N;
con_func = @(theta)Target.stable_con(theta);
for k = 1 : multi_startpt_number
    while true
        theta = randn(1, nvar);
        if  con_func(theta) < 0
            break;
        end
    end
    theta_set(k, :) = theta;
end
% nosie set
sim_N = 1000;
mean_N = 100;

ddd = randn(sim_N, 2, mean_N);
% ddd = [sim_N, 2, mean_N];

%% simple
smiple_opt = struct();
smiple_opt.FDM_width = 1e-8;
smiple_opt.FDM_shrink_size = 0.5; % less than 1.0
smiple_opt.FDM_inf_width = 1e-12;
smiple_opt.STEP_width = 1e-6;
smiple_opt.STEP_shrink_size = 0.5; % less than 1.0
smiple_opt.STEP_inf_width = 1e-12;
smiple_opt.STEP_alpha = 0.8;% less than 1.0
smiple_opt.MaxIteration = 100;

% ResultChk = simple_gradient(Target, theta_set, ddd, opt);
% [final_value, idx] = find_min(ResultChk.value);

%%
inter_opt = struct();
inter_opt.FDM_width_f = 1e-8;
inter_opt.FDM_width_c = 1e-8;
inter_opt.MaxIteration = 100;
inter_opt.mu = 1;% barriar parameter
inter_opt.alpha = 1e-6;
inter_opt.alpha_shrink = 0.5;


%%
stable_flag = @(theta) Target.stable_con(theta) < 0;
stable_con = @(theta) Target.stable_con(theta);
eval_func = @(theta, ddd) Target.eval_func(theta, ddd);

ResultChk = struct();

IterationChk = 0;

theta = theta_set(1, :);


% A single constraint
Lambda = eye(1);% dual variable
% BFGS initial Hessian
B_f = eye(length(theta));
B_c = eye(1);


FDM_grad_f = nan(size(theta));
FDM_grad_c = nan(size(theta));
FDM_width_f = inter_opt.FDM_width_f;
FDM_width_c = inter_opt.FDM_width_c;
parfor itr1 = 1 : length(theta)
    one_hot = zeros(size(theta));
    one_hot(itr1) = 1;
    tilde = FDM_width_f;
    upper_tilde_theta = theta + tilde*one_hot;
    lower_tilde_theta = theta - tilde*one_hot;
    FDM_grad_f(itr1) = (eval_func(upper_tilde_theta,ddd)-eval_func(lower_tilde_theta,ddd))/(2*tilde);
    fprintf('Target function %dth FDM grad complete !!\n', itr1);
end
parfor itr1 = 1 : length(theta)
    one_hot = zeros(size(theta));
    one_hot(itr1) = 1;
    tilde = FDM_width_c;
    upper_tilde_theta = theta + tilde*one_hot;
    lower_tilde_theta = theta - tilde*one_hot;
    FDM_grad_c(itr1) = (stable_con(upper_tilde_theta)-stable_con(lower_tilde_theta))/(2*tilde*stable_con(theta));
    fprintf('Constraint function %dth FDM grad complete !!\n', itr1);
end
lambda = -inter_opt.mu/stable_con(theta);
% inter update
A = [B_f, -FDM_grad_c'; lambda*FDM_grad_c, B_c];
b = [-FDM_grad_f'+-FDM_grad_c'*lambda; inter_opt.mu-B_c*lambda];
p = A\b;
idx =  0;
disp('start update')
while true
    delta = inter_opt.alpha*(inter_opt.alpha_shrink^(idx));
    new_theta = theta + delta*p(1:length(theta))';
    if stable_flag(new_theta) 
        theta = new_theta;
        lambda = lambda + delta*p(length(theta)+1:end)';
        disp('end update')
        break;
    end
    if delta < sqrt(eps)
        disp('Impossible update')
        break;
    end
    idx = idx + 1;
end
%% evaluation
% % x_multi = ResultChk(idx).theta;
% % 
% % rng(data.eval_seed)
% % evaluation_N = 100;
% % ddd = randn(sim_N, 2, evaluation_N);
% % f_slqr_set = zeros(evaluation_N, 1);
% % f_elqr_set = zeros(evaluation_N, 1);
% % f_multi_set = zeros(evaluation_N, 1);
% % yyy_slqr_set = zeros(sim_N, 3, evaluation_N);
% % yyy_elqr_set = zeros(sim_N, 3, evaluation_N);
% % yyy_multi_set = zeros(sim_N, 3, evaluation_N);
% % 
% % sim_t = (0:sim_N-1)'*data.Ts;
% % Ts = data.Ts;
% % sys_slqr = sys_all_slqr({'y_node1','u_controlled1'}, {'d_node1'});
% % sys_elqr = sys_all_elqr({'y_node1','u_controlled1'}, {'d_node1'});
% % evaluation_func = @(theta, noise)Target.eval_func(theta, noise);
% % parfor k = 1 : evaluation_N
% %     yyy_slqr_set(:,:,k) = lsim(sys_slqr, ddd(:,:,k), sim_t, 'zoh');
% %     f_slqr_set(k) = norm(yyy_slqr_set(:,:,k));
% %     yyy_elqr_set(:,:,k) = lsim(sys_elqr, ddd(:,:,k), sim_t, 'zoh');
% %     f_elqr_set(k) = norm(yyy_elqr_set(:,:,k));
% %     [f_multi_set(k), ~, yyy_multi_set(:,:,k)] = evaluation_func(x_multi, ddd(:,:,k));
% % end
% % 
% % disp('reward')
% % mean(f_multi_set)
% % mean(f_elqr_set)
% % mean(f_slqr_set)
% % 
% % sim_t = (0:sim_N-1)'*Target.Ts;
% % figure,
% % plot(sim_t, mean(yyy_multi_set(:,2,:),3), 'g');
% % hold on
% % plot(sim_t, mean(yyy_elqr_set(:,2,:),3), 'r');
% % plot(sim_t, mean(yyy_slqr_set(:,2,:),3), 'b');

%% local
function [v,i] = find_min(varargin)
    [v,i] = min(cell2mat(varargin));
end

function ResultChk = simple_gradient(Target, theta_set, ddd, opt)
    stable_flag = @(theta) Target.stable_con(theta) < 0;
    eval_func = @(theta, ddd) Target.eval_func(theta, ddd);

    ResultChk = struct();
    for multi_point = 1 : size(theta_set, 1)

        IterationChk = 0;
        theta = theta_set(multi_point, :);
        FDM_width = opt.FDM_width;
        FDM_shrink_size = opt.FDM_shrink_size;
        FDM_inf_width = opt.FDM_inf_width;
        while true
            FDM_grad = nan(size(theta));
            parfor itr1 = 1 : length(theta)
                one_hot = zeros(size(theta));
                one_hot(itr1) = 1;
                idx = 0;
                while true
                    tilde = FDM_width*(FDM_shrink_size^(idx));
                    if tilde < FDM_inf_width
                       tilde = FDM_inf_width; 
                    end
                    upper_tilde_theta = theta + tilde*one_hot;
                    lower_tilde_theta = theta - tilde*one_hot;
                    if stable_flag(upper_tilde_theta) && stable_flag(lower_tilde_theta) 
                        FDM_grad(itr1) = (eval_func(upper_tilde_theta,ddd)-eval_func(lower_tilde_theta,ddd))/(2*tilde);
                        break;
                    end
                    if tilde == FDM_inf_width
                       break; 
                    end
                    idx = idx + 1;
                end
            end

            if any(isnan(FDM_grad))
                disp('FDM Impossile !!');
                break;
            end

            idx = 0;
            while true 
                STEP = opt.STEP_width*(opt.STEP_shrink_size^(idx));
                if STEP < opt.STEP_inf_width
                   STEP = opt.STEP_inf_width; 
                end
                new_theta = STEP*FDM_grad + theta;
                if stable_flag(new_theta)
                    theta = opt.STEP_alpha*STEP*FDM_grad + theta;
                    break; 
                end
                if STEP == opt.STEP_inf_width
                    break;
                end
                idx = idx + 1;
            end
            if STEP == opt.STEP_inf_width
                disp('Reach STEP infinium width')
            end
            if IterationChk == opt.MaxIteration
                break;
            end
            IterationChk = IterationChk + 1;
        end
        eval_value = eval_func(theta, ddd);
        ResultChk(multi_point).value = eval_value;
        ResultChk(multi_point).theta = theta;
        fprintf('%dth parameter : %d iteration update\n', multi_point, IterationChk);
    end
    disp('End MultiStart Local Optimization');
end
