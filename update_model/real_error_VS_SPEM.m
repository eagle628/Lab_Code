clearvars -except Store
Store = load('C:\Users\NaoyaInoue\Desktop\figure_set\node30_confirm_ss_oe_spem\node1\np0\data');

%% Extract Requaired Valriavle
sys_local = Store.sys_local;
sys_env = Store.sys_env;
[a,b,c,d] = ssdata(sys_env);

number = 58;
id_dim = order(Store.sys_ID_set{number});

K2 = [-1,-1];
K3 = [-1];

[sys_spem_spem1, sys_real_spem1] = G_checkz_d(sys_env, sys_local, Store.sys_ID_set{number}, K2, K3);
[sys_spem_oe1, sys_real_oe1] = G_checkz_d(sys_env, sys_local, Store.sys_oe_set{number}, K2, K3);
[sys_spem_oe_ss1, sys_real_oe_ss1] = G_checkz_d(sys_env, sys_local, Store.sys_oe_ss_set{number}, K2, K3);
[sys_spem_bal1, sys_real_bal1] = G_checkz_d(sys_env, sys_local, balred(sys_env,id_dim), K2, K3);


[sys_spem_spem2, sys_real_spem2] = G_checkz_d(sys_env, sys_local, Store.sys_ID_set{number});
[sys_spem_oe2, sys_real_oe2] = G_checkz_d(sys_env, sys_local, Store.sys_oe_set{number});
[sys_spem_oe_ss2, sys_real_oe_ss2] = G_checkz_d(sys_env, sys_local, Store.sys_oe_ss_set{number});
[sys_spem_bal2, sys_real_bal2] = G_checkz_d(sys_env, sys_local, balred(sys_env,id_dim));

% K2 = [-1,1];
% K3 = [0];
% [sys_spem_spem2, sys_real_spem2] = compare_func(sys_env, sys_local, Store.sys_ID_set{number}, K2, K3);
% [sys_spem_oe2, sys_real_oe2] = compare_func(sys_env, sys_local, Store.sys_oe_set{number}, K2, K3);
% [sys_spem_oe_ss2, sys_real_oe_ss2] = compare_func(sys_env, sys_local, Store.sys_oe_ss_set{number}, K2, K3);
% [sys_spem_bal2, sys_real_bal2] = compare_func(sys_env, sys_local, balred(sys_env,id_dim), K2, K3);
% 
% 
% K2 = [-1,-1];
% K3 = [-1];
% [sys_spem_spem3, sys_real_spem3] = compare_func(sys_env, sys_local, Store.sys_ID_set{number}, K2, K3);
% [sys_spem_oe3, sys_real_oe3] = compare_func(sys_env, sys_local, Store.sys_oe_set{number}, K2, K3);
% [sys_spem_oe_ss3, sys_real_oe_ss3] = compare_func(sys_env, sys_local, Store.sys_oe_ss_set{number}, K2, K3);
% [sys_spem_bal3, sys_real_bal3] = compare_func(sys_env, sys_local, balred(sys_env,id_dim), K2, K3);

%% local function
function [sys_spem, sys_real] = G_checkz_d(sys_env, sys_local, sys_id, K2, K3)
if sys_id.Ts ~= 0
    sys_id = ss(d2c(sys_id));
end
if nargin < 4
    Q = diag([1, 1000]);
    R = 1e-3;
    [~, K, ~, ~] = Retrofit.design(sys_local, sys_id, Q, R);
    K2 = K(1:2);
    K3 = K(3:end)*sys_id.B;
end

Ae = sys_env.A;
Be = sys_env.B;
Ce = sys_env.C;
De = sys_env.D;

A = sys_local.A;
L = sys_local({'y'},{'v'}).B;
B = sys_local({'y'},{'u'}).B;
R = sys_local({'y'},{'u'}).B;
% R = [0;1];
C = sys_local({'y'},{'u'}).C;
Gamma = sys_local({'w'},{'u'}).C;
S = sys_local({'y'},{'v'});
[a,b,c,d] = ssdata(S);
S = ss(a,b,c(2,:),d(2,:));

Ae_apx = sys_id.A;
Be_apx = sys_id.B;
Ce_apx = sys_id.C;
De_apx = sys_id.D;


% % % % % % % w_hat & w & vを生成するシステム
l_dim = order(sys_local);
id_dim = order(sys_id);
e_dim = order(sys_env);

A_xi = [
        A+L*De*Gamma, L*Ce, L*(De-De_apx)*Gamma, -L*Ce_apx;
        Be*Gamma, Ae, Be*Gamma, zeros(e_dim,id_dim);
        zeros(l_dim,l_dim), zeros(l_dim,e_dim), A+L*De_apx*Gamma+B*K2*C+B*K3*Gamma, L*Ce_apx;
        zeros(id_dim,l_dim), zeros(id_dim,e_dim), Be_apx*Gamma, Ae_apx;
        ];
B_xi = [
        zeros(l_dim,1);
        zeros(e_dim,1);
        R;
        zeros(id_dim,1);
        ];
C_xi = [
        Gamma, zeros(1,e_dim), Gamma, zeros(1,id_dim);
        De*Gamma, Ce, De*Gamma, zeros(1,id_dim);
        zeros(1,l_dim+e_dim), Gamma, zeros(1,id_dim);
        ];
    
% G_w_v_wh__d = minreal(ss(A_xi, B_xi, C_xi, []));
G_w_v_wh__d = ss(A_xi, B_xi, C_xi, []);
G_w_v_wh__d.OutputGroup.('w') = 1;
G_w_v_wh__d.OutputGroup.('v') = 2;
G_w_v_wh__d.OutputGroup.('w_hat') = 3;

% % % % % % 予期せぬフィードバック項
G_feed_apx = ss();
G_feed_apx.A = [
                A+L*De_apx*Gamma, L*Ce_apx;
                Be_apx*Gamma, Ae_apx;
                ];
G_feed_apx.B = [
                L;
                zeros(id_dim,1);
                ];
G_feed_apx.C = [De_apx*Gamma, Ce_apx];
G_feed_apx.D = eye(1);

G_feed_true = ss();
G_feed_true.A = [
                A+L*De*Gamma, L*Ce;
                Be*Gamma, Ae;
                ];
G_feed_true.B = [
                L;
                zeros(e_dim,1);
                ];
G_feed_true.C = [De*Gamma, Ce];
G_feed_true.D = eye(1);

sys_spem = minreal(S*G_feed_apx*(G_w_v_wh__d('v')-sys_id*G_w_v_wh__d('w')));
sys_real = minreal(S*G_feed_true*(sys_env-sys_id)*G_w_v_wh__d('w_hat'));

% sys_spem = S*G_feed_apx*(G_w_v_wh__d('v',:)-sys_id*G_w_v_wh__d('w',:));
% sys_real = S*G_feed_true*(sys_env-sys_id)*G_w_v_wh__d('w_hat',:);


end

function sys = sys_series(sys1, sys2)
    [a1, b1, c1, d1] = ssdata(sys1);
    [a2, b2, c2, d2] = ssdata(sys2);
    
    a = [
        a1, tool_zeros(a1,a2);
        b2*c1, a2;
        ];
    b = [
        b1;
        b2*d1;
        ];
    c = [
        d2*c1, c2;
        ];
    d = d2*d1;
    
    sys = ss(a, b, c, d);
end

function Zero_mat = tool_zeros(M1, M2)
    Zero_mat = zeros(size(M1,1), size(M2,2));
end
