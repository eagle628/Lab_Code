clear
close all

%%

%% Exsample GA (https://jp.mathworks.com/help/gads/ga.html)
% xi = linspace(-6,2,300);
% % % yi = linspace(-4,4,300);
% % [X,Y] = meshgrid(xi,yi);
% % Z = ps_example([X(:),Y(:)]);
% % Z = reshape(Z,size(X));
% % surf(X,Y,Z,'MeshStyle','none')
% % colormap 'jet'
% % view(-26,43)
% % xlabel('x(1)')
% % ylabel('x(2)')
% % title('ps\_example(x)')
% % 
% % %%
% % opt = optimoptions('ga','PlotFcn',@gaplotbestf);
% % 
% % rng default % For reproducibility
% % [x, fval] = ga(@ps_example,2, [], [], [], [], [], [], [], opt)
%% 
data.node_number = 2;
data.net_seed = 6;
data.net_local = 1;
data.Ts = 0.1;

data.controller_n = 4;
%% define model
% % % network version
net = network_swing_simple(data.node_number, [1,2], [2,12]*1e-3, 1, [0.1,6], 0.8, data.net_seed);
model = swing_network_model(net, data.net_local, data.Ts);
% SLQR
net.add_controller(data.net_local, diag([1,1]), 1);
sys_all_slqr = net.get_sys();
sys_all_slqr = net.add_io(sys_all_slqr, data.net_local, 'node1');
sys_all_slqr = net.get_sys_controlled(sys_all_slqr);
net.controllers = {};
% ELQR
[~, sys_env] = net.get_sys_local(data.net_local);
net.add_controller(data.net_local, sys_env, diag([1,1]), 1);
sys_all_elqr = net.get_sys();
sys_all_elqr = net.add_io(sys_all_elqr, data.net_local, 'node1');
sys_all_elqr = net.get_sys_controlled(sys_all_elqr);
net.controllers = {};

%% define controller model
ss_model = gen_ss_tridiag(data.controller_n, model.ny, model.nu);
%
% env_n = 2;
% env_model = gen_ss_tridiag(env_n, 1, 1);
% ss_model = gen_ss_Eretro_controller(env_model, c2d(model.sys_local, model.Ts, model.discrete_type));

%% pre
rng default % For reproducibility
sim_N = 1000;
mean_N = 100;
% ddd = randn(sim_N, 2, mean_N);
ddd = [sim_N, 2, mean_N];

env = model.RL_env_all_prime;
AE = env.A;
BE = env(:, 'u_node1').B;
CE = env({'yhat','what'},:).C;
R  = env(:, 'd_node1').B;
S  = env('y', :).C;

local = model.sys_local_discrete;
[Ap,Bp,Cp,~] = ssdata(local);

%% GA
% opt = optimoptions(...
%                     'ga',...
%                     'UseParallel',true,...
%                     'MaxStallGenerations',1e3,...
%                     'Display', 'iter',...
%                     'PopulationSize',100 ...
%                     );
% 
% 
% [x, fval] = ga(@(theta)eval_func(theta, AE,BE,CE,R,S,data.Ts, ss_model, ddd), ss_model.N, [], [], [], [], [], [], @(theta)stable_con(theta, Ap,Bp,Cp, ss_model), opt)


%% PSO
multi_startpt_number = 1000;
theta_set = zeros(multi_startpt_number, ss_model.N);

nvar = ss_model.N;
parfor k = 1 : multi_startpt_number
    while true
        theta = randn(1, nvar);
        if  stable_con(theta, Ap,Bp,Cp, ss_model) < 0
            break;
        end
    end
    theta_set(k, :) = theta;
end
tpoints = CustomStartPointSet(theta_set);


problem = createOptimProblem('fmincon',...
    'objective',@(theta)eval_func(theta, AE,BE,CE,R,S,data.Ts, ss_model, ddd),...
    'nonlcon',@(theta)stable_con(theta, Ap,Bp,Cp, ss_model),...
    'x0',theta_set(1, :), ...
    'options',...
    optimoptions(...
                @fmincon,...
                'Algorithm','interior-point',...
                'Display','iter'...
                )...
    );

fmincon(problem);
ms = MultiStart('UseParallel',true,'Display','iter');
[x,fval,eflag,output,manymins] = run(ms, problem, tpoints);
%% evaluation
rng(6)
evaluation_N = 100;
ddd = randn(sim_N, 2, evaluation_N);
f_slqr_set = zeros(evaluation_N, 1);
f_elqr_set = zeros(evaluation_N, 1);
f_global_set = zeros(evaluation_N, 1);
yyy_slqr_set = zeros(sim_N, 3, evaluation_N);
yyy_elqr_set = zeros(sim_N, 3, evaluation_N);
yyy_global_set = zeros(sim_N, 3, evaluation_N);

sim_t = (0:sim_N-1)'*data.Ts;
Ts = data.Ts;
sys_slqr = sys_all_slqr({'y_node1','u_controlled1'}, {'d_node1'});
sys_elqr = sys_all_elqr({'y_node1','u_controlled1'}, {'d_node1'});
parfor k = 1 : evaluation_N
    yyy_slqr_set(:,:,k) = lsim(sys_slqr, ddd(:,:,k), sim_t, 'zoh');
    f_slqr_set(k) = norm(yyy_slqr_set(:,:,k));
    yyy_elqr_set(:,:,k) = lsim(sys_elqr, ddd(:,:,k), sim_t, 'zoh');
    f_elqr_set(k) = norm(yyy_elqr_set(:,:,k));
    [f_global_set(k), ~, yyy_global_set(:,:,k)] = eval_func(x, AE,BE,CE,R,S,Ts, ss_model, ddd(:,:,k));
end
%% local
function [f, grad, y] = eval_func(theta, AE,BE,CE,R,S,Ts, controller, ddd)
    controller.set_params(theta);
    [Ak,Bk,Ck,Dk,dAk,dBk,dCk,dDk] = controller.get_ss();

    Anew = [AE+BE*Dk*CE, BE*Ck; Bk*CE, Ak];
    Bnew = [R; tools.zeros(Ak, R)];
    Cnew = [S, tools.zeros(S, Ak); Dk*CE, Ck];
    sys = ss(Anew, Bnew, Cnew, [], Ts);

    f = 0;
    grad = zeros(length(theta), 1);
    if length(ddd) == 3
        ddd = randn(ddd);
    end
    for iter1 = 1 : size(ddd, 3)
        [y,~,x] = lsim(sys, ddd(:,:,iter1), []);
        f = f + norm(y);
        if nargin == 2 
            for iter2 = 1 : length(dAk)
                dAnew = [AE+BE*dDk{iter2}*CE, BE*dCk{iter2}; dBk{iter2}*CE, dAk{iter2}];
                dBnew = [R; tools.zeros(Ak, R)];
                dCnew = [S, tools.zeros(S, Ak); dDk{iter2}*CE, dCk{iter2}];
                dsys = ss(Anew, [dAnew,dBnew], Cnew, [dCnew,tools.zeros(dCnew, dBnew)], Ts);
                dy = lsim(dsys, [x, ddd(:,:,iter1)], []);
                grad(iter2) = grad(iter2) + sum(sum(y.*dy));
            end
        end
    end
    f = f/size(ddd, 3);
end

function [c, ceq] = stable_con(theta, Ap,Bp,Cp, controller)
    controller.set_params(theta);
    [Ak,Bk,Ck,Dk] = controller.get_ss();
    A_all = [Ak,Bk*Cp; Bp*Ck, Ap+Bp*Dk*Cp];
    pole = eig(A_all);
    c= max(abs(pole))- 1;
    ceq = [];
end

function out = ReLU(x)
    if x < 0
        out = 0;
    else
        out = x^2;
    end
end
