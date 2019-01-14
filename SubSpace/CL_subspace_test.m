% section 5.5.3
% Four Step CLIVC4 algorithm
%% main
clear
close all
%% Prepare Plan and Controller and Noise Filter
rng(20)
GGG = 3;
while true
    G = tf([rand(1,GGG-1),0],[1,rand(1,GGG-1)]);
    G = rss(GGG-1);
    G.D = randn(1);
    if isstable(G)
        break;
    end
end

na = GGG-1;
nb = GGG;

% na = 5;
% nb = 6;
% 
iii = 1;
% while true
rng(iii)
CCC = 2;
while true
    C = tf([rand(1,1)],[1,rand(1,CCC-2),0]);
    feed = loopsens(G, C);
    if feed.Stable
        break;
    end
%     iii = iii + 1;
end
% C = sys_local_vw;

rng(1024)
HHH = GGG-1;
while true
    H = tf([1,randn(1,HHH)],[1,randn(1,HHH)]);
    if isstable(H)
        break;
    end
end
% H = ss(G);
% H.B = zeros(GGG-1, 1);
% H.B(4) = 1;

% H_d = c2d(H, Ts, 'foh');

nc = HHH;
nd = HHH;
% nc = 5;
% nd = 5;


%% Generate Response : Availabel r1 r2 y u
disp('Generate I/O data')
% sample config
N = 10000;
Ts = 0.01;
% disvretiozation
G_d = c2d(G, Ts, 'foh');
C_d = c2d(C, Ts, 'foh');
H_d = c2d(H, Ts, 'foh');
% genrate mesurementa data
r1 = randn(N, 1)*0;
r2 = randn(N, 1)*1;
d  = randn(N, 1)*0;
r  = r1 + lsim(C_d, r2);
cloop_d = loopsens(G_d,C_d);
yyy = lsim(G_d*cloop_d.Si, r);
ddd = lsim(cloop_d.So*H_d, d);
y = yyy + ddd;
u = r - lsim(C_d, y);

%% main

sys = moesp_DT(iddata(zeros(N,1),[u,y], Ts), 3, 30);

s = 30;
n = 2;

[ydat, ny] = size(y);
[udat, nu] = size(u);

N = N - 2*s + 1;
R_1  = blkhankel(r, 1, N, s);
R_s1 = blkhankel(r, s+1, N, s);
U_1  = blkhankel(u, 1, N, s);
U_s1 = blkhankel(u, s+1, N, s);
Y_s1 = blkhankel(y, s+1, N, s);

% LQ decomposition
L = triu(qr([R_1; R_s1; U_1; U_s1; Y_s1]'))';
L = L(1:s*(4*nu+ny),:);

L51 = L((nu*s)*4+1:end, 1:nu*s);
L52 = L((nu*s)*4+1:end, nu*s+1:2*nu*s);

L31 = L((nu*s)*2+1:(nu*s)*3, 1:nu*s);
L32 = L((nu*s)*2+1:(nu*s)*3, nu*s+1:2*nu*s);
L41 = L((nu*s)*3+1:(nu*s)*4, 1:nu*s);
L42 = L((nu*s)*3+1:(nu*s)*4, nu*s+1:2*nu*s);

P1 = L31 - (L31*L41'+L32*L42')*(L41*L41'+L42*L42')*L41;
P2 = L32 - (L31*L41'+L32*L42')*(L41*L41'+L42*L42')*L42;

Gh = (L51*P1'+L52*P2')/(P1*P1'+P2*P2')^(1/2);

% SVD 
[U1, S1, V1] = svd(Gh);
% sigular value
singular = diag(S1);

% C and A
U1_L = U1(:, 1:n);%*diag(sqrt(singular(1:n)));
C = U1_L(1:ny, :);
A = U1_L(1:ny*(s-1), :)\U1_L(ny+1:s*ny, :);

% B & D
U1_R = U1(:, n+1:end);%*diag(sqrt(singular(n+1:end)));
U1_R = U1_R';

m = size(U1_R,1);
alpha = zeros(s*m, s);
for itr = 0 : s-1
    alpha(m*itr+1:m*(itr+1), 1:end-itr) = U1_R(:, 1+itr:end);
end
beta = U1_R*(L51*L41'+L52*L42')/(L41*L41'+L42*L42');
beta = reshape(beta, numel(beta), 1);

DB = blkdiag(eye(ny),U1(1:ny*(s-1), :)) \ (alpha \ beta);

D = DB(1:ny, :);
B = DB(ny+1:ny+n, :);

idsys = ss(A, B, C, D, Ts);

%% local



function sys = CLTS4SID(data, n, s)
    u = data.u;
    y = data.y;
    r = data.r;
    Ts= dat.Ts;
    N = data.N;

    N = N - 2*s + 1;

    [ydat, ny] = size(y);
    [udat, nu] = size(u);
    [rdat, nr] = size(r);

    U_s1 = blkhankel(u, s+1, N, s);
    R_s1 = blkhankel(r, s+1, N, s);
    R_1  = blkhankel(r, 1, N, s);

    Z = ort_pro(U_s1, R_s1, R_1);

    Z_1_Ns  = Z(:, 1:end-s);
    Z_s1_Ns = Z(:, s+1:end);
    Y_s1_Ns = blkhankel(y, s+1, N-s, s);

    % LQ decomposition
    R = triu(qr([Z_1_Ns;Z_s1_Ns;Y_s1_Ns]'))';
    R = R(1:s*(ny+nu), :);

    % SVD
    R22 = R(s*nu+1:end, s*nu+1:end);
    [U1, S1, ~]=svd(R22);

    % sigular value
    singular = diag(S1);

    % C and A
    Ok = U1(:, 1:n)*diag(sqrt(singular(1:n)));
    C=Ok(1:ny, :);
    A=Ok(1:ny*(s-1), :)\Ok(ny+1:s*ny, :);

    % B and D
    L1 = U1(:, n+1:end)';
    R11 = R(1:s*nu, 1:s*nu);
    R21 = R(s*nu+1:end, 1:s*nu);
    M1 = L1*R21/R11;
    m = ny*s-n;
    M = zeros(m*s, nu);
    L = zeros(m*s, ny+n);
    for k=1:s
        M((k-1)*m+1:k*m, :)=M1(:, (k-1)*nu+1:k*nu);
        L((k-1)*m+1:k*m, :)=[L1(:, (k-1)*ny+1:k*ny) L1(:, k*ny+1:end)*Ok(1:end-k*ny,:)];
    end
    DB=L\M;
    D=DB(1:ny, :);
    B=DB(ny+1:end, :);

    sys = ss(A, B, C, D, Ts);
end


%% local local
function H = blkhankel(u, i, j, s)
    nu = size(u, 2);
    H = zeros(s*nu, j);
    for itr = 1 : nu
        H(1:nu:end, :) = hankel(u(i:i+s-1, itr), u(i+s-1:i+j+s-2, itr ));
    end
end

function PI = Pi1(A)
    PI = A'*((A*A')\A);
end

function PI = Pi2(A)
    PI = eye(size(A,2)) - A'*((A*A')\A);
end

% Orthognal projection
function Z = ort_pro(A, B, C)
    Z = A * ( Pi1(B) + Pi2(B)*C'*((C*Pi2(B)*C')\(C*Pi2(B))) );
end
