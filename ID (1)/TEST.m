close all
clear

P_s = tf([9,0],[1,1,9]);
H_s = tf(1,[1,8]);

Ts = 0.1;
P = c2d(P_s,Ts);
H = c2d(H_s,Ts);

N = 100000;
r = (rand(N,1)-0)*1;
% d = (rand(N,1)-0.5)*1;
cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',1);
d = cn3();
clear cn3

% u = lsim(P,r) + lsim(G_ud,d);
u = r;
y = lsim(P,r) + lsim(H,d);

idmodel = pem_discrete(y,u,Ts);
init = c2d(rss(6),Ts);

mine = idmodel.pem([6,7,6,6,6,0],init);
mine_tf = mine.B/mine.A/mine.F;

% mine = idmodel.armax([6,7,6,0]);
% mine_tf = mine.B/mine.A;

% mine = idmodel.oe([7,6,0]);
% mine_tf = mine.B/mine.F;

model_pem = pem(iddata(y,u,Ts),idss(init),'Display','on');

bode(mine_tf)
hold on
bode(model_pem)
bode(P)

legend('my','kumikomi','Plant')