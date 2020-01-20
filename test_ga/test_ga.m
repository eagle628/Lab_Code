clear
% close all

%% condition
data.node_number = 2;
data.net_seed = 8;
data.net_local = 1;
data.Ts = 0.1;

data.controller_n = 4;

data.eval_seed = 6;

data.ga_seed = 12;
data.pso_seed = 12;
data.train_input_seed = 28;
data.controller_struct_enable = false;
data.train_input_form = true;
% data.controller_gain_enable = true;

%% For test
% controller_struct_enables = [false, true];
% train_input_forms = [false, true];
% train_input_seeds = [1,2,3,4,5];
% opt_seeds = [6,7,8];

% for controller_struct_enable = controller_struct_enables
%     for train_input_form = train_input_forms
%         for opt_seed = opt_seeds
%             if train_input_form
%                 data.ga_seed = opt_seed;
%                 data.pso_seed = opt_seed;
%                 data.controller_struct_enable = controller_struct_enable;
%                 data.train_input_form = train_input_form;
%                 data.input_seed = 0;
%                 test_func(data);
%             else
%                 for train_input_seed = train_input_seeds
%                     data.ga_seed = opt_seed;
%                     data.pso_seed = opt_seed;
%                     data.controller_struct_enable = controller_struct_enable;
%                     data.train_input_form = train_input_form;
%                     data.input_seed = train_input_seed;
%                     test_func(data);
%                 end
%             end
%         end
%     end
% end

% function test_func(data)
%% define model
% % % network version
net = network_swing_simple(data.node_number, [1,2], [2,12]*1e-3, 1, [0.1,6], 0.8, data.net_seed);
net.Adj_ref = net.Adj_ref*0;
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
%% pre
rng(data.train_input_seed)
sim_N = 1000;
mean_N = 10;
if data.train_input_form
    ddd = randn(sim_N, 2, mean_N);
else
    ddd = [sim_N, 2, mean_N];
end
%% GA
opt = optimoptions(...
                    'ga',...
                    'UseParallel',true,...
                    'MaxStallGenerations',1e3,...
                    'Display', 'iter',...
                    'PopulationSize',100 ...
                    );

rng(data.ga_seed) % For reproducibility
% [x_ga, fval_ga] = ga(@(theta)Target.eval_func(theta, ddd), ss_model.N, [], [], [], [], [], [], @(theta)Target.stable_con(theta), opt);
%% fmincon multi start
multi_startpt_number = 100;
theta_set = zeros(multi_startpt_number, ss_model.N);


rng(data.pso_seed)
nvar = ss_model.N;
con_func = @(theta)Target.stable_con(theta);
parfor k = 1 : multi_startpt_number
    while true
        theta = randn(1, nvar);
        if  con_func(theta) < 0
            break;
        end
    end
    theta_set(k, :) = theta;
end
tpoints = CustomStartPointSet(theta_set);


problem = createOptimProblem('fmincon',...
    'objective',@(theta)Target.eval_func(theta, ddd),...
    'nonlcon',@(theta)Target.stable_con(theta),...
    'x0',theta_set(1, :), ...
    'options',...
    optimoptions(...
                @fmincon,...
                'Algorithm','interior-point',...
                'Display','iter' ...
                )...
    );

x1 = fmincon(problem);

problem = createOptimProblem('fmincon',...
    'objective',@(theta)Target.eval_func(theta, ddd),...
    'nonlcon',@(theta)Target.stable_con(theta),...
    'x0',theta_set(1, :), ...
    'options',...
    optimoptions(...
                @fmincon,...
                'Algorithm','interior-point',...
                'Display','iter', ...
                'SpecifyObjectiveGradient', true, ...
                'SpecifyConstraintGradient', true ...
                )...
    );
x2 = fmincon(problem);
% ms = MultiStart('UseParallel',true,'Display','iter');
% [x_multi,fval_multi,eflag_multi,output_multi,manymins_multi] = run(ms, problem, tpoints);
%% PSO
% rng(1024)
% opt = optimoptions(...
%                 'particleswarm', ...
%                 'Display','iter', ...
%                 'FunctionTolerance', 1e-6, ...
%                 'UseParallel',true, ...
%                 'SwarmSize', 1000 ...
%                 );
%
% target = @(theta)Target.PSOcon(...
%                 Target.eval_func(theta, ddd), ...
%                 Target.stable_con(theta) ...
%                 );
% [x_pso,fval_pso,exitflag_pso,output_pso] = particleswarm(target, ss_model.N, [], [], opt);
%% evaluation
% tmp = canon(balreal(c2d(balred(balreal(sys_env),2),data.Ts,'zoh')),'modal');
% tmp_theta = diag(tmp.A);
% tmp_theta = [tmp_theta;diag(tmp.A,-1)];
% tmp_theta = [tmp_theta;diag(tmp.A,1)];
% tmp_theta = [tmp_theta;tmp.B];
% tmp_theta = [tmp_theta;tmp.C'];
% tmp_theta = [tmp_theta;tmp.D];
% tmp_theta = tmp_theta';
% x_multi = tmp_theta;
% x_ga = tmp_theta;
% load(sprintf('test_rl_grad_belief_%d_epi1.mat',belief_N))
% x_ga = x_';
% x_multi = x_';
% load 2019-12-07=23-07-54 x_ga x_pso
% x_multi = x_pso;

sim_N = 500;
rng(347)
evaluation_N = 100;
ddd = zeros(sim_N, 2, evaluation_N);
tmp  = load('2020-01-10=00-15-41');
sigma = tmp.train.opt_policy.target.pi_sigma;
parfor itr1 = 1: evaluation_N
    rng(itr1)
%     ddd(:,2,itr1) = (sigma^2)*randn(sim_N, 1, 1);
    ddd(:,:,itr1) = randn(sim_N, 2, 1);
end
f_slqr_set = zeros(evaluation_N, 1);
f_elqr_set = zeros(evaluation_N, 1);
f_multi_set = zeros(evaluation_N, 1);
f_ga_set = zeros(evaluation_N, 1);
yyy_slqr_set = zeros(sim_N, 3, evaluation_N);
yyy_elqr_set = zeros(sim_N, 3, evaluation_N);
yyy_multi_set = zeros(sim_N, 3, evaluation_N);
yyy_ga_set = zeros(sim_N, 3, evaluation_N);

sim_t = (0:sim_N-1)'*data.Ts;
Ts = data.Ts;
sys_slqr = sys_all_slqr({'y_node1','u_controlled1'}, {'d_node1'});
sys_elqr = sys_all_elqr({'y_node1','u_controlled1'}, {'d_node1'});
evaluation_func = @(theta, noise)Target.eval_func(theta, noise);
parfor_progress(evaluation_N);
parfor k = 1 : evaluation_N
    yyy_slqr_set(:,:,k) = lsim(sys_slqr, ddd(:,:,k), sim_t, 'zoh');
    f_slqr_set(k) = sum(sum(yyy_slqr_set(:,:,k).^2));
    yyy_elqr_set(:,:,k) = lsim(sys_elqr, ddd(:,:,k), sim_t, 'zoh');
    f_elqr_set(k) = sum(sum(yyy_elqr_set(:,:,k).^2));
    [f_multi_set(k), ~, yyy_multi_set(:,:,k)] = evaluation_func(x_multi, ddd(:,:,k));
    f_multi_set(k) = sum(sum(yyy_multi_set(:,:,k).^2));
    [f_ga_set(k), ~, yyy_ga_set(:,:,k)] = evaluation_func(x_ga, ddd(:,:,k));
    f_ga_set(k) = sum(sum(yyy_ga_set(:,:,k).^2));
    parfor_progress();
end
parfor_progress(0);

disp('reward')
mean(f_multi_set)
mean(f_elqr_set)
mean(f_slqr_set)
mean(f_ga_set)

% sim_t = (0:sim_N-1)'*Target.Ts;
% figure,
% plot(sim_t, mean(yyy_multi_set(:,2,:),3), 'g');
% hold on
% plot(sim_t, mean(yyy_elqr_set(:,2,:),3), 'r');
% plot(sim_t, mean(yyy_slqr_set(:,2,:),3), 'b');
% plot(sim_t, mean(yyy_ga_set(:,2,:),3), 'k');
% 
% figure,
% boxplot(...
% [f_ga_set,f_multi_set,f_slqr_set,f_elqr_set],...
% {'GA','Multi','SLQR','ELQR'}...
% )

% savetime = char(datetime);
% savetime = strrep(savetime,'/','-');
% savetime = strrep(savetime,':','-');
% savetime = strrep(savetime,' ','=');
% save(savetime)

% end
%%

% while true
% x = randn(1, Target.controller.N);
%     if Target.stable_con(x) < 0
%         break;
%     end
% end
%
%%
episode = 10;
% load(sprintf('data_rl_n%d_grad_b%d/test_rl_grad_belief_%d_epi%d.mat',model.net.N,belief_N,belief_N, episode))
% load(sprintf('data_rl_n%d_grad_b%d/test_rl_grad_n%d_belief_%d_epi%d.mat',model.net.N,belief_N,model.net.N,belief_N, episode))
load(sprintf('test_rl_grad_belief_%d_epi%d.mat',belief_N, episode))
% load(sprintf('test_rl_grad_n%d_belief_%d_epi%d.mat',model.net.N,belief_N, episode))
x_test = x_';
input_ = [input_',zeros(size(input_'))];

step_sizes = [1e-8, -1e-8];
opt_eval = Target.eval_func(x_test, input_);
area_grad = zeros(1, ss_model.N, length(step_sizes));
for ivar = 1 : ss_model.N
one_hot = zeros(1, ss_model.N);
one_hot(ivar) = 1;
    for k = 1 : length(step_sizes)
        area_grad(1, ivar, k) = Target.eval_func(x_test+one_hot*step_sizes(k), input_)-opt_eval;
    end
end
cent_grad = (area_grad(:,:,1)-area_grad(:,:,2))./2e-8;
cent_grad = cent_grad./norm(cent_grad)

figure
plot(grad(1:end-1)./norm(grad(1:end-1)))
hold on
plot(-cent_grad)
legend('RL','FDM')
title(sprintf('Episode-%d',episode))
%% 
close all
sss1 = load('tmp1');
sss2 = load('tmp2');

figure('Name','Controller_0_Multi')
boxplot(...
[sss2.f_multi_set,sss1.f_multi_set,sss1.f_slqr_set,sss1.f_elqr_set],...
{'Multi1','Multi2','SLQR','ELQR'})

xlabel('制御器設計方法')
ylabel('評価値')
title('制御器設計手法による評価値の分布')
f = gcf;
f.Units = 'centimeters';
f.Position = [0,0,15,10.5];
ylim([0 inf]);
f.Children.YGrid = 'on';
f.Children.FontSize = 16;
f.Children.FontName = 'SanSelif';
f.Children.YMinorGrid = 'on';

figure('Name','Controller_0_GA')
boxplot(...
[sss2.f_ga_set,sss1.f_ga_set,sss1.f_slqr_set,sss1.f_elqr_set],...
{'GA1','GA2','SLQR','ELQR'})
xlabel('制御器設計方法')
ylabel('評価値')
title('制御器設計手法による評価値の分布')
f = gcf;
f.Units = 'centimeters';
f.Position = [0,0,15,10.5];
ylim([0 inf]);
f.Children.YGrid = 'on';
f.Children.FontSize = 16;
f.Children.FontName = 'SanSelif';
f.Children.YMinorGrid = 'on';

sss1 = load('tmp3');
sss2 = load('tmp4');

figure('Name','Controller_1_Multi')
boxplot(...
[sss2.f_multi_set,sss1.f_multi_set,sss1.f_slqr_set,sss1.f_elqr_set],...
{'Multi1','Multi2','SLQR','ELQR'})

xlabel('制御器設計方法')
ylabel('評価値')
title('制御器設計手法による評価値の分布')
f = gcf;
f.Units = 'centimeters';
f.Position = [0,0,15,10.5];
ylim([0 inf]);
f.Children.YGrid = 'on';
f.Children.FontSize = 16;
f.Children.FontName = 'SanSelif';
f.Children.YMinorGrid = 'on';

figure('Name','Controller_1_GA')
boxplot(...
[sss2.f_ga_set,sss1.f_ga_set,sss1.f_slqr_set,sss1.f_elqr_set],...
{'GA1','GA2','SLQR','ELQR'})
xlabel('制御器設計方法')
ylabel('評価値')
title('制御器設計手法による評価値の分布')
f = gcf;
f.Units = 'centimeters';
f.Position = [0,0,15,10.5];
ylim([0 inf]);
f.Children.YGrid = 'on';
f.Children.FontSize = 16;
f.Children.FontName = 'SanSelif';
f.Children.YMinorGrid = 'on';

