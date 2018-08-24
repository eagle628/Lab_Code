clear
% Rnadma Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
n = network_swing_simple(30, [1, 2], [2, 10]*1e-2, 1, [1, 5], [], 10);

n.Adj_ref = n.Adj_ref*0;
Q = diag([1, 1000]);
R = 1e-3;


sys_org = n.get_sys_controlled({2}, {'n2'});
    % n.add_io(1),n.add_io(3)をやったのと同値
N = 50000;
Ts = 0.01;
d = randn(N, 2);
d2 = randn(N, 2);
t = (0:N-1)'*Ts;
% 外乱dに対するvとｗの実際（シミュレーション）の応答
v = lsim(sys_org({'v_n2'}, 'd_n2'), d, t);
w = lsim(sys_org({'w_n2'}, 'd_n2'), d, t);

sys_all = n.get_sys();
%環境モデルの抽出
[sys_local, sys_env] = n.get_sys_local(1);

n_model = 6;
m = model_cont_std(n_model, 0);%何次モデルにするか
sys_lower = balred(sys_env, 6);
m.set_tf(sys_lower);% 低次元化パラメータを初期値にセット
m.lsim_type = 'foh';
m.fit_mymqt(t, w, v);% マルカール法によるパラメータ推定
model1 = minreal(ss(m));% 最小実現


n.add_controller(1, Q, R);
sys_c1 = n.get_sys_controlled({1, 3}, {'n1', 'n3'});

n.controllers = {};
n.add_controller(1, model1, Q, R);

sys_c2 = n.get_sys_controlled({1, 3}, {'n1', 'n3'});

close all
% Insert impulse Noise for Original System
[y0, t0] = impulse(sys_org({'y_n1'}, 'd_n1'), 0:0.01:2000);
% Insert impulse Noise for retro fit controller of Node1
[y1, t1] =impulse(sys_c1({'y_n1', 'xhat_c1', 'u_c1'}, 'd_n1'), 0:0.01:2000);
% Insert impluse Noise for Extend retro fit controller of Node1
[y2, t2] = impulse(sys_c2({'y_n1', 'xhat_c1', 'u_c1'}, 'd_n1'), 0:0.01:2000);
h=tools.mysubplot(t0, {y0(:, 2, 1), nan(size(t0)), nan(size(t0))}, [600, 600], 'k', 'LineWidth', 1.5);
tools.myplot2(h(1),t1, y1(:, 2, 1), 'b', 'LineWidth', 1.5)
tools.myplot2(h(1),t2, y2(:, 2, 1), 'r', 'LineWidth', 1.5)
tools.myplot2(h(3),t1, y1(:, 4, 1), 'b', 'LineWidth', 1.5)
tools.myplot2(h(3),t2, y2(:, 4, 1), 'r', 'LineWidth', 1.5)
tools.myplot2(h(2),t1, y1(:, 2, 1)-y1(:, 4, 1), 'b', 'LineWidth', 1.5)
tools.myplot2(h(2),t2, y2(:, 2, 1)-y2(:, 4, 1), 'r', 'LineWidth', 1.5)
xlim([0, 500])


[y0_a, t0_a] = impulse(sys_org({'y_n3'}, 'd_n1'), 0:0.01:2000);
[y1_a, t1_a] =impulse(sys_c1({'y_n3'}, 'd_n1'), 0:0.01:2000);
[y2_a, t2_a] = impulse(sys_c2({'y_n3'}, 'd_n1'), 0:0.01:2000);

[fig, h]=tools.myplot([600, 600], t0_a, y0_a(:, 2, 1), 'k', 'LineWidth', 1.5);
tools.myplot2(h,t1_a, y1_a(:, 2, 1), 'b', 'LineWidth', 1.5)
tools.myplot2(h,t2_a, y2_a(:, 2, 1), 'r', 'LineWidth', 1.5)