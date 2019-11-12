close all
clear

sys_order = 6;
K_order = 2;


Ts = 0.1;
N = 10;
t = (0:N-1)'*Ts;
u = randn(N, 1);

%% Open loop
% % x0 = zeros(sys_order, 1);
% % sys = rss(sys_order);
% % [y1,t1,x1] = lsim(sys, u, t, x0, 'zoh');
% % sys_d = c2d(sys, Ts, 'zoh');
% % [y2,t2,x2] = lsim(sys_d, u, t, x0);
% % 
% % 
% % [a, b, c, d] = ssdata(sys_d);
% % abcd = [a, b;c, d];
% % y3 = zeros(N, 1);
% % x3 = zeros(size(x2));
% % x3(1, :) = x0;
% % for k = 1: N
% %     xy = abcd*[x3(k, :)'; u(k, :)'];
% %     x3(k+1, :) = xy(1:sys_order);
% %     y3(k, :) = xy(sys_order+1:end);
% % end

%% Closed loop
dis_type = 'zoh';
sys = rss(sys_order);
sys.d = 1;
P  = c2d(sys, Ts, dis_type);

rng(10)
iter = 1;
while true
    K = ss(randn(K_order), randn(K_order, 1), randn(1, K_order), randn(1), Ts);
    loop = feedback(P, K, +1);
    if isstable(loop)
        break;
    end
    disp(iter)
    iter = iter + 1;
end
    
[y1,t,x1] = lsim(loop, u, t, zeros(order(loop),1));

x_p = zeros(sys_order, N);
x_k = zeros(K_order, N);
y_p = zeros(1, N);
y_k = zeros(1, N);

[ap,bp,cp,~] = ssdata(P);
[ak,bk,ck,dk] = ssdata(K);
% 
loop2 = ss([ap,bp*ck;bk*cp,ak], [bp;zeros(K_order,1)], [cp,zeros(1,K_order)], [], Ts);
[y2,t,x2] = lsim(loop2, u, t, zeros(order(loop2),1),'zoh');

% 
for k = 1 : N
    y_p(:, k) = cp*x_p(:, k);
    y_k(:, k) = ck*x_k(:, k) + dk*y_p(:, k);
    if k ~= N
        x_k(:, k+1) = ak*x_k(:, k) + bk*y_p(:, k);
        x_p(:, k+1) = ap*x_p(:, k) + bp*(y_k(:, k)+u(k, :));
    end
end
