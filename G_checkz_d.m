%% NOTE
% 引数が4の時：Kが与えらえられた状態でかつそのKは状態フィードバックのもの
% 引数が5の時：K2,K3はそれぞれローカルシステムの次数と環境システムの入出力次元数に当たる
% 引数が3の時：sys_idを環境システムとして，勝手にLQRでゲインが確定する．

%%
function [sys_spem, sys_real] = G_checkz_d(sys_env, sys_local, sys_id, K)
    if sys_id.Ts ~= 0
        sys_id = ss(d2c(sys_id));
    end
    [a,b,c,d] = ssdata(sys_id);
    sys_id = ss(a,b,c,d);
    if nargin < 4
        Q = diag([1, 1]);
        R = 1e-3;
        [~, K, ~, ~] = Retrofit.design(sys_local, sys_id, Q, R);
    end
    K2 = K(1:2);

    Ae = sys_env.A;
    Be = sys_env.B;
    Ce = sys_env.C;
    De = sys_env.D;

    A = sys_local.A;
    L = sys_local({'y'},{'v'}).B;
    B = sys_local({'y'},{'u'}).B;
    R = sys_local({'y'},{'u'}).B;% omegaにノイズが導入されたとき
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

%     new_local_sys = Retrofit.generate_fb(sys_local, sys_id);
%     new_local_sys = new_local_sys('y','u');
%     K3 = K(3:end)*new_local_sys.B;
%     K3 = K(3:end)*sys_id.B;
    K3 = sys_id.C*K(3:end)';

    % % % % % % % w_hat & w & vを生成するシステム
    l_dim = order(sys_local);
    id_dim = order(sys_id);
    e_dim = order(sys_env);

    A_xi = [
            A+L*De*Gamma, L*Ce, L*(De-De_apx)*Gamma, -L*Ce_apx;
            Be*Gamma, Ae, Be*Gamma, zeros(e_dim,id_dim);
            zeros(l_dim,l_dim), zeros(l_dim,e_dim), A+L*De_apx*Gamma+B*K2*C+B*K3*Gamma, L*Ce_apx;
%             zeros(l_dim,l_dim), zeros(l_dim,e_dim), A+L*De_apx*Gamma+tmp, L*Ce_apx;
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

    sys_spem = minreal(S*G_feed_apx*(G_w_v_wh__d('v')-sys_id*G_w_v_wh__d('w')), [], false);
    sys_real = minreal(S*G_feed_true*(sys_env-sys_id)*G_w_v_wh__d('w_hat'), [], false);

    % sys_spem = S*G_feed_apx*(G_w_v_wh__d('v',:)-sys_id*G_w_v_wh__d('w',:));
    % sys_real = S*G_feed_true*(sys_env-sys_id)*G_w_v_wh__d('w_hat',:);
end


%% 
%{
XR = [
eye(2, 2), sys_local({'y'},{'v'})*feedback(iv_IDsys_set{1}.G,-sys_local({'w'},{'v'})), -sys_local({'y'},{'v'})*feedback(1, -iv_IDsys_set{1}.G*sys_local({'w'},{'v'}));
zeros(1, 2), feedback(1,-sys_local({'w'},{'v'})*iv_IDsys_set{1}.G), feedback(-sys_local({'w'},{'v'}),iv_IDsys_set{1}.G)
];
%}
