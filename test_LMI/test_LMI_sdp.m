clear
close all

m = 1;
n = 6;

F1 = [];

xi_lower = sdpvar(m, n, 'full');
xi_upper = sdpvar(m, n, 'full');

c_mean = (xi_lower+xi_upper)./2;
c_val  = xi_upper - c_mean;

lambda = sdpvar(m, n, 'full');
F1 = [F1, vec(lambda)>=0];

M1 = diag(sum(lambda.*(c_val.^2 - c_mean.^2), 1));
tmp = -lambda.*c_mean;
U = [];
for itr1 = 1 : m
    U = [U, diag(tmp(itr1, :))];
end
M4 = diag(vec(-lambda'));

M = [M1,U;U',M4];

W = kron(eye(m), ones(1, n));

sys = rss(n,1,m);

A = sys.A;
B = sys.B;

P = sdpvar(n, n, 'symmetric','complex');
F1 = [F1, P>=0];

gamma = 1;

O1 = [A'*P+P*A, P*B*W; W'*B'*P, zeros(m*n)];
O2 = blkdiag(eye(n)./gamma, zeros(m*n));
O3 = M;

O = O1 + O2 + O3;
S = [P*B; zeros(m*n, m)];

SDP = [O, S; S', -gamma*eye(m)];

F1 = [F1, SDP<0];

sol1 = optimize(F1, [])


F2 = [];
x = randn(n,1);
q = sdpvar(m, n, 'full');
x = repmat(x', m, 1);

phi = (double(c_val).^2 - double(c_mean).^2).*(x.^2) + 2*double(c_mean).*q.*x - q.^2;

F2 = [F2, vec(phi)>0 ];

sol2 = optimize(F2, [])

disp('xi_lower')
disp(double(xi_lower))
disp('xi_upper')
disp(double(xi_upper))

q_ans = vec(double(q'));

w = W*q_ans
