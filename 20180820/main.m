close all
clear

P_s = tf([9,0],[1,1,9]);
K_s = tf(1,[1,8]);

Ts = 0.1;
P = c2d(P_s,Ts);
K = c2d(K_s,Ts);

S = 1/(1+P*K);
G_yr = P*S;
G_ur = S;

G_yd = S;
G_ud = -K*S;

N = 10000;
r = rand(N,1)*1;
d = rand(N,1)*1;

u = lsim(G_ur,r) + lsim(G_ud,d);
y = lsim(G_yr,r) + lsim(G_yd,d);

data = iddata(y,u,Ts);
opt = armaxOptions;
% d_ID = armax(data,[2,3,2,0]);
d_ID = arx(data,[2,3,0]);

n = 2;
nk = 1;
a = sym('a', [n, 1]);
b = sym('b', [n+1, 1]);

ak = sym('ak', [nk, 1]);
bk = sym('bk', [nk+1, 1]);


syms q

Phat = q.^(n:-1:0)*b/(q.^(n:-1:0)*[1;a])
Khat = q.^(nk:-1:0)*bk/(q.^(nk:-1:0)*[1;ak])

Shat = 1/(1+Phat*Khat)

[numK, denK] = tfdata(K, 'v');
numK = numK/denK(1);
denK = denK/denK(1);

KK = subs(Khat, [ak; bk], [denK(2:end), numK]')
[num, den] = numden(KK);
num = double(coeffs(num, q))
den = double(coeffs(den, q))
KKK = tf(num, den, [])


% [num, den] = numden(simplify(S))
