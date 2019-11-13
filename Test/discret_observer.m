close all
clear

n = 6;
p = 3;
m = 2;
sys = drss(n,p,m);

[A,B,C,~] = ssdata(sys);
Ts = 0.1;
sys = ss(A,B,C,[],Ts);
Q = eye(n);
R = eye(p);
L = -dlqr(A', C', Q, R)';

N = 1000;
t = (0:N-1)'*Ts;

u = randn(m, N);
x_true = nan(n, N);
x_hat = nan(n, N);
y = nan(p, N);

x_true(:, 1) = zeros(n, 1);
x_hat(:, 1) = zeros(n, 1);


Ahat = A;
Bhat = B;
Chat = C;

for k = 1:N
    y(:, k) = C*x_true(:, k);
    x_true(:, k+1) = A*x_true(:, k) + B*u(:, k);
    x_hat(:, k+1) = (Ahat+L*Chat)*x_hat(:, k) + [Bhat, -L]*[u(:,k);y(:,k)];
end
