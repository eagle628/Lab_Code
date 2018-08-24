%{
lsim は任意入力に対する応答を計算できる，
odeは自立システムの厳密解（数値解）を計算する．
つまり，odeは外部入力を受け付けず，システム内部で明示的に入力を示せる場合，
（多くの場合には，状態フィードバック）についてしか応答を計算でじゅばう，

lsimは，外部入力を許容する分システムを離散化して，厳密解では，ないくしているといえるかもしれない．


%}

close all
clear

G = tf(1,[1,1]);
G = ss(G);

N = 100;
Ts = 0.001;
t_s = 0:Ts:Ts*(N-1);
t_s = t_s;

u = zeros(size(t_s));
y_l = lsim(G, u, t_s, G.B);

[ y_i, t_i] = impulse(G,t_s(end));

figure
plot(t_s,y_l,'-*');
hold on
plot(t_i,y_i,'-*');

legend('lsim','impulse')

 