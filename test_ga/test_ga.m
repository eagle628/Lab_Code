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
% ss_model = gen_ss_tridiag(data.controller_n, model.ny, model.nu);
%
env_n = 2;
env_model = gen_ss_tridiag(env_n, 1, 1);
ss_model = gen_ss_Eretro_controller(env_model, c2d(model.sys_local, model.Ts, model.discrete_type));


%%
opt = optimoptions('ga',...
                    'UseParallel',true,...
                    'MaxStallGenerations',1e3,...
                    'Display', 'iter',...
                    'PopulationSize',100 ...
                    );

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


[x, fval] = ga(@(theta)eval_func(theta, AE,BE,CE,R,S,data.Ts, ss_model, ddd), ss_model.N, [], [], [], [], [], [], @(theta)stable_con(theta, Ap,Bp,Cp, ss_model), opt)

%% evaluation
ddd = randn(sim_N, 2);
yyy_slqr = lsim(sys_all_slqr({'y_node1','u_controlled1'}, {'d_node1'}), ddd, (0:sim_N-1)'*data.Ts, 'zoh');
f_slqr = norm(yyy_slqr)
yyy_elqr = lsim(sys_all_elqr({'y_node1','u_controlled1'}, {'d_node1'}), ddd, (0:sim_N-1)'*data.Ts, 'zoh');
f_elqr = norm(yyy_elqr)
[f_global, yyy_global] = eval_func(x, AE,BE,CE,R,S,data.Ts, ss_model, ddd)
%% local
function [f, y] = eval_func(theta, AE,BE,CE,R,S,Ts, controller, ddd)
    controller.set_params(theta);
    [Ak,Bk,Ck,Dk] = controller.get_ss();

    Anew = [AE+BE*Dk*CE, BE*Ck; Bk*CE, Ak];
    Bnew = [R; tools.zeros(Ak, R)];
    Cnew = [S, tools.zeros(S, Ak); Dk*CE, Ck];
    sys = ss(Anew, Bnew, Cnew, [], Ts);

    f = 0;
    if length(ddd) == 3
        ddd = randn(ddd);
    end
    for k = 1 : size(ddd, 3)
        y = lsim(sys, ddd(:,:,k), []);
        f = f + norm(y);
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
