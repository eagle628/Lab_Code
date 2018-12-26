clear
close all
%%
eig(28)
test_dim = 4;
G = rss(test_dim);
H = tf(1);

N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;

u = randn(N,1);
d = randn(N,1)*0;

yyy = lsim(G, u, t, [], 'foh');
ddd = lsim(H, d, t, [], 'foh');

y = yyy + ddd;

G_d = c2d(G, Ts, 'foh');

%%
n = test_dim;% identification dim

[ndat,ny]=size(y);
[mdat,nu]=size(u);
if ndat ~= mdat
    error('Y and U have different length.')
end

s = 30;
% block Hankel matrix
N = ndat-s+1;
Y = zeros(s*ny,N);
U = zeros(s*nu,N);
%     sN = sqrt(N);
%     sy = y'/sN;
%     su = u'/sN;
sy = y';
su = u';
for s=1:s
    Y((s-1)*ny+1:s*ny, :)=sy(:, s:s+N-1);
    U((s-1)*nu+1:s*nu, :)=su(:, s:s+N-1);
end

% LQ decomposition
R=triu(qr([U;Y]'))';
R=R(1:s*(ny+nu), :);

% SVD
R22 = R(s*nu+1:end, s*nu+1:end);
[U1, S1, ~]=svd(R22);

% sigular value
singular = diag(S1);

% C and A
Ok = U1(:, 1:n)*diag(sqrt(singular(1:n)));
% C=Ok(1:ny, :);
% A=Ok(1:ny*(d-1), :)\Ok(ny+1:d*ny, :); 

%%
Q = sdpvar(n,n,'full');
P = sdpvar(n,n,'symmetric');
% Q = semivar(n,n,'full');
% P = semivar(n,n,'symmetric');

Objective = norm(Ok(1:ny*(s-1), :)*Q-Ok(ny+1:s*ny, :)*P,'fro');

ep = 1e-6;
Constraints = [];
Constraints = [Constraints, P>ep];
alpha = 1*eye(2);
beta = [0,1;0,0];
M_func = @(P,Q,alpha,beta) kron(alpha,P) + kron(beta,Q) + kron(beta',Q');
Constraints = [Constraints, M_func(P,Q,alpha,beta)>=ep];
Constraints = [Constraints, trace(P) == n];

% ops = sdpsettings('solver','cplex','verbose',0);

sol = optimize(Constraints,Objective);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
 solution_P = value(P);
 solution_Q = value(Q);
else
 disp('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end

A = solution_Q\solution_P
eig(A)
G_d.A
eig(G_d.A)
