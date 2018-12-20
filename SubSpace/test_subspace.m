clear
close all

G = rss(2);
H = tf(1);

N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;

u = randn(N,1);
d = randn(N,1)*0;

yyy = lsim(G, u, t, [], 'foh');
ddd = lsim(H, d, t, [], 'foh');

y = yyy + ddd;

%%
idsys_DT = moesp_DT(iddata(y, u, Ts), 2, 30);
idsys_CT = moesp_CT(iddata(y, u, Ts), 2, 30);

bode(G,idsys_DT,idsys_CT)
legend('original','MOESP-DT','MOESP-CT','location','best')




