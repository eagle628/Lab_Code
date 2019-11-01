clear
close all

seed1 = 8;
Node_number = 2;
net = network_swing_simple(Node_number, [1,2], [2,12]*1e-3, 1, [0.1,6], 0.8, seed1);
% net.Adj_ref = net.Adj_ref*0;
% local system
c_n = 1;
% set model
Ts = 0.1;
model = swing_network_model(net, c_n, Ts);

Q = diag([1,1]);
R = 1;

sys_all = model.sys_all;

ob_y = {'y_node1'};
in_d = {'d_node1'};

%% 
Te = 10;
ini = [0,0];

%% original (model current dynamics)
MPC_mode = false;
[ori_local_x_all1, ori_env_x_all1, ori_rect_x_all1, ori_mpc_u_all1, ori_t1, ori_d_L1] = test_sim(model, Te, ini, Q,R,MPC_mode);

ori_sys2 = sys_all(ob_y, in_d);
initial = zeros(order(ori_sys2), 1);
initial(1:2) = ini;
ori_local_x_all2 = lsim(ori_sys2, ori_d_L1, ori_t1, initial,'zoh');
ori_local_x_all3 = lsim(ori_sys2, ori_d_L1, ori_t1, initial,'foh');
ori_local_x_all4 = lsim(c2d(ori_sys2,model.Ts,'zoh'), ori_d_L1, ori_t1, initial);
ori_local_x_all5 = lsim(c2d(ori_sys2,model.Ts,'foh'), ori_d_L1, ori_t1, initial);

figure('Name','Original')
plot(ori_local_x_all1,'r')
hold on
plot(ori_local_x_all2,'b')
plot(ori_local_x_all3,'g')
plot(ori_local_x_all4,'c')
plot(ori_local_x_all5,'k')

%% lqr (model current dynamics)
MPC_mode = true;
[lqr_local_x_all1, lqr_env_x_all1, lqr_rect_x_all1, lqr_mpc_u_all1, lqr_t1, lqr_d_L1] = test_sim(model, Te, ini, Q,R,MPC_mode);

net.add_controller(model.c_n, Q, R);
lqr_sys2 = net.get_sys_controlled(sys_all);
lqr_sys2 = lqr_sys2(ob_y, in_d);
net.controllers = {};

initial = zeros(order(lqr_sys2), 1);
initial(1:2) = ini;
lqr_local_x_all2 = lsim(lqr_sys2, lqr_d_L1, lqr_t1, initial,'zoh');
lqr_local_x_all3 = lsim(lqr_sys2, lqr_d_L1, lqr_t1, initial,'foh');
lqr_local_x_all4 = lsim(c2d(lqr_sys2,model.Ts,'zoh'), lqr_d_L1, lqr_t1, initial);
lqr_local_x_all5 = lsim(c2d(lqr_sys2,model.Ts,'foh'), lqr_d_L1, lqr_t1, initial);

figure('Name','LQR')
plot(lqr_local_x_all1,'r')
hold on
plot(lqr_local_x_all2,'b')
plot(lqr_local_x_all3,'g')
plot(lqr_local_x_all4,'c')
plot(lqr_local_x_all5,'k')

%% autensity input () 
% u = randn(length(lqr_d_L1),1);
% [local_x_all, env_x_all, rect_x_all, mpc_u_all,t,d_L] = test_sim2(model, Te, ini, u);
% test_sys = sys_all(ob_y, {'u_node1'});
% 
% xxx = lsim(test_sys, u, t, [ini,0,0],'zoh');
%% 
rect_sys = Retrofit.generate_rectifier(model.sys_local, ss(0));
sys_all = model.sys_all({model.port_y,model.port_w,model.port_v}, {model.port_control,model.port_d_L});
% RL_env = rect_sys*sys_all;

[A2,B2,C2,D2] = ssdata(rect_sys);
[A1, B1, C1, D1] = ssdata(sys_all);
Ak = [A1, tools.zeros(A1, A2); B2*C1, A2];
Bk = [B1; B2*D1];
Ck = [D2*C1 C2];
Dk = D2*D1;
RL_env_all = ss(Ak,Bk,Ck,Dk);
RL_env_all.OutputGroup = rect_sys.OutputGroup;
RL_env_all.InputGroup = sys_all.InputGroup;


RL_env = RL_env_all({'yhat','what'},:);
lqr_sys = RL_env_all({'x','y'},:);
lqr_sys = c2d(lqr_sys, Ts, 'foh');
RL_env = c2d(RL_env, Ts, 'foh');
[test_x_all1,y_w_all, mpc_u_all, t, dL] = new_dynamics(model, lqr_sys, RL_env, Te, ini, Q, R, true);

net.add_controller(model.c_n, Q, R);
lqr_sys3 = net.get_sys_controlled(model.sys_all);
lqr_sys3 = lqr_sys3(ob_y, in_d);
net.controllers = {};

initial = zeros(order(lqr_sys3), 1);
initial(1:2) = ini;
test_x_all2 = lsim(lqr_sys3, dL, t, initial,'zoh');
test_x_all3 = lsim(lqr_sys3, dL, t, initial,'foh');
test_x_all4 = lsim(c2d(lqr_sys3,model.Ts,'zoh'), dL, t, initial);
test_x_all5 = lsim(c2d(lqr_sys3,model.Ts,'foh'), dL, t, initial);
figure('Name','New dynamics')
plot(test_x_all1(:,1:2),'r');
hold on
plot(test_x_all2,'b')
plot(test_x_all3,'g')
plot(test_x_all4,'c')
plot(test_x_all5,'k')

%% new dynamics feed connect
% belief_N(What numbers do observe signal store ?)
belief_N = 30;
%  define init controller
seed_define = 2;
rng(seed_define)
controller_n = 4;% state
controller_l = 1;% output
nnn = model.ny+model.nw;
controller_m = belief_N*nnn;% input
A1 = diag(ones(1,(belief_N-1)*nnn),-nnn);
B1 = zeros(size(A1,1),nnn);
B1(1:nnn,1:nnn) = eye(nnn);
C1 = A1;
D1 = B1;
recorder = ss(A1,B1,C1,D1,model.Ts);
recorder.InputGroup.yhat = [1,2];
recorder.InputGroup.what = 3;
P = model.sys_local_discrete({'y','w'},{'u'});
iter = 1;
while true
    controller = drss(controller_n, controller_l, controller_m);
    loop = feedback(P, controller*recorder, +1);
    if isstable(loop)
        break;
    end
    disp(iter)
    iter = iter+1;
end
[A2,B2,C2,D2] = ssdata(controller);
Ak = [A1, tools.zeros(A1, A2); B2*C1, A2];
Bk = [B1; B2*D1];
Ck = [D2*C1 C2];
Dk = D2*D1;
true_controller = ss(Ak,Bk,Ck,Dk,Ts);
true_controller.InputGroup = recorder.InputGroup;
true_controller.OutputGroup.u = 1;

test_sys = feedback(c2d(RL_env_all({'yhat','what','y'},{'u_node1','d_node1'}),Ts,'zoh'), true_controller, 1,[1,2,3],+1);
%% local
function [x_all,y_w_all, mpc_u_all, t, d_L] = new_dynamics(model, lqr_sys, RL_env, Te, ini, Q, R, MPC_mode)
    sim_N = Te/model.Ts + 1;
    t = (0:model.Ts:Te)';
    d_L = randn(sim_N, 2);
    % get memory
    x_all = zeros(sim_N, model.local_nx+model.env_nx+model.rect_nx);
    y_w_all = zeros(sim_N, model.ny+model.nw);
    mpc_u_all = zeros(sim_N, model.nu);
    % set initial
    x_all(1, 1:2) = ini';
    % calculate MBC gain
    K = dlqr(model.A, model.B, Q, R);
    if ~MPC_mode
        K = zeros(size(K));
    end
    K1 = K(1:model.local_nx);
    K2 = K(model.local_nx+1 : end);
    K = [K1,-K2,-K1];
    for k = 1 : sim_N-1
        % LQR Controller 
        xy = lqr_sys.C*x_all(k, :)';
        mpc_u_all(k, :) = (K*xy)';
        %  observe S_(k+1)
        x_all(k+1, :)= (RL_env.A*x_all(k, :)' + RL_env.B*[mpc_u_all(k, :)'; d_L(k, :)'])';
        y_w_all(k, :) = (RL_env.C*x_all(k, :)')';
    end 
end



function [local_x_all, env_x_all, rect_x_all, mpc_u_all,t,d_L] = test_sim(model, Te, ini, Q,R,MPC_mode)
    sim_N = Te/model.Ts + 1;
    t = (0:model.Ts:Te)';
    d_L = randn(sim_N, 2);
    local_x_all = zeros(sim_N, model.local_nx);
    env_x_all = zeros(sim_N, model.env_nx);
    rect_x_all = zeros(sim_N, model.rect_nx);
    y_w_v_all = zeros(sim_N, model.ny+model.nw+model.nv);
    mpc_u_all = zeros(sim_N, model.nu);
    % set initial
    local_x_all(1, :) = ini';
    % calculate MBC gain
    K = dlqr(model.A, model.B, Q, R);
    if ~MPC_mode
        K = zeros(size(K));
    end
        K1 = K(1:model.local_nx);
        K2 = K(model.local_nx+1 : end);
    for k = 1 : sim_N-1
         % LQR Controller 
        ywv = model.dynamics(local_x_all(k, :)', env_x_all(k, :)');
        y_w_v_all(k, :) = ywv';
        mpc_u_all(k, :) = ([K1,-K2]*rect_x_all(k, :)')' -(K1*y_w_v_all(k, 1:model.ny)')';
        %  observe S_(k+1)
        [~, local_ne_x, env_ne_x] = model.dynamics(local_x_all(k, :)', env_x_all(k, :)', mpc_u_all(k, :), d_L(k, :)');
        local_x_all(k+1, :) = local_ne_x';
        env_x_all(k+1, :) = env_ne_x';
%         [rect_ne_x, rect_yw] = model.rect_dynamics(local_x_all(k, :)', env_x_all(k, :)', rect_x_all(k, :)');
        [rect_ne_x, rect_yw] = local_rect(model, local_x_all(k, :)', env_x_all(k, :)', rect_x_all(k, :)',  ywv);
        rect_x_all(k+1, :) = rect_ne_x';
    end 
end

 function [ne_rect_x, rect_yw] = local_rect(model, local_x, env_x, rect_x, ywv)
    all_x = model.L_and_E_to_all(local_x, env_x);
%     ywv = model.observe_all*all_x;
    rect_xy = model.system_rect*[rect_x;ywv];
    ne_rect_x = rect_xy(1:model.rect_nx, :);
    rect_yw = rect_xy(model.rect_nx+1:end, :);
 end

 function [local_x_all, env_x_all, rect_x_all, mpc_u_all,t,d_L] = test_sim2(model, Te, ini, u)
    sim_N = Te/model.Ts + 1;
    t = (0:model.Ts:Te)';
    d_L = zeros(sim_N, 2);
    local_x_all = zeros(sim_N, model.local_nx);
    env_x_all = zeros(sim_N, model.env_nx);
    rect_x_all = zeros(sim_N, model.rect_nx);
    y_w_v_all = zeros(sim_N, model.ny+model.nw+model.nv);
    mpc_u_all = zeros(sim_N, model.nu);
    % set initial
    local_x_all(1, :) = ini';
    for k = 1 : sim_N-1
         % LQR Controller 
        ywv = model.dynamics(local_x_all(k, :)', env_x_all(k, :)');
        y_w_v_all(k, :) = ywv';
        mpc_u_all(k, :) = u(k)
        %  observe S_(k+1)
        [~, local_ne_x, env_ne_x] = model.dynamics(local_x_all(k, :)', env_x_all(k, :)', mpc_u_all(k, :), d_L(k, :)');
        local_x_all(k+1, :) = local_ne_x';
        env_x_all(k+1, :) = env_ne_x';
%         [rect_ne_x, rect_yw] = model.rect_dynamics(local_x_all(k, :)', env_x_all(k, :)', rect_x_all(k, :)');
        [rect_ne_x, rect_yw] = local_rect(model, local_x_all(k, :)', env_x_all(k, :)', rect_x_all(k, :)',  ywv);
        rect_x_all(k+1, :) = rect_ne_x';
    end 
end
