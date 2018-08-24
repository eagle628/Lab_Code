clear 
close all
s = tf('s');

Ps = 1/(s+2);
Cs = (s+8)*(s+9)/((s+3)*(s+5)*(s+6)*(s+7));
%Cs = 1/(s+3);
Hs = 1;
Vs = 1;

Ts = 0.1;
P = c2d(Ps,Ts);
C = c2d(Cs,Ts);
H = Hs;
V = Vs;

G_eu = C*H/(1-C*P);
G_ey = H/(1-C*P);

G_ru = V/(1-C*P);
G_ry = P*V/(1-C*P);

N = 5000;
rng(10);
e = rand(N,1)*1;
r = rand(N,1)*0;
t = (0:N-1)'*Ts;

u = lsim(G_eu,e,t)+lsim(G_ru,r,t);
y = lsim(G_ey,e,t)+lsim(G_ry,r,t);

data = iddata(y,u,Ts);

nk = 1;
%model_oe = oe(data,[2,2,0]);
model_armax = armax(data,[1,1,5,nk]);

bode(P,inv(C),model_armax)
legend('P','C_inv','identification')