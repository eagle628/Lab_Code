function [D_Sys2] = my_BBM_SISO(u,y,system,Ts)
%% Explain Function 
% Black Box Model 
% Single Output SingleInput
% Delay とりあえず無視
%% 
na = system(1);
nb = system(2);
nc = system(3);
nd = system(4);
nf = system(5);
nk = system(6);
%
N = length(u);
% 初期条件を考慮
d = nk-1;
Y = zeros(number,na);
U = zeros(number,nb);
W = zeros(number,nf);
E = zeros(number,nc);
V = zeros(number,nd);

for i = 1 : na
    Y(:,i) = [zeros(i,1);y(1:end-i)];
end
for i = 1 : nb
    U(:,i) = [zeros(i+d,1);u(1:end-d-i)];
end


PHI = [-Y,U];
theata2 = PHI\y;

D_Sys2 = tf([zeros(1,nk),theata2(na+1:na+nb)'],[1,theata2(1:na)'],Ts,'Variable','z^-1');
end
