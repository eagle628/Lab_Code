function [D_Sys2] = myARX(u,y,system,Ts)

na = system(1);
nb = system(2);
nk = system(3);
%{
number = length(u);
H = zeros((na+nb)*(number-na),na+nb);
phi = zeros((na+nb)*(number-na),na+nb);
delta_y = zeros((na+nb)*(number-na),1);
for i = 1+max(na,nb)+nk:number
    phi_y = zeros(na,1);
    phi_u = zeros(nb,1);
    for k = 1:na
        phi_y(k) = -y(i-k);
    end
    for k = 1:nb
        phi_u(k) = u(i-k-nk+1);
    end
    phi_l = [phi_y;phi_u;];
    phi((i-1-na)+1,:) = phi_l';
    H((na+nb)*(i-1-na)+1:(na+nb)*(i-na),:) = phi_l*phi_l';
    delta_y((na+nb)*(i-1-na)+1:(na+nb)*(i-na)) = phi_l'*y(i);
end

theata1 = (H'*H)\(H'*delta_y);
D_Sys1 = tf([zeros(1,nk),theata1(na+1:na+nb)'],[1,theata1(1:na)'],Ts,'Variable','z^-1');
%}

%
number = length(u);
% 初期条件を考慮
d = nk-1;
Y = zeros(number,na);
U = zeros(number,nb);
for i = 1 : na
    Y(:,i) = [zeros(i,1);y(1:end-i)];
end
for i = 1 : nb
    U(:,i) = [zeros(i+d,1);u(1:end-d-i)];
end

PHI = [-Y(na+1:end,:),U(na+1:end,:)];
theata2 = PHI\y(na+1:end);

D_Sys2 = tf([zeros(1,nk),theata2(na+1:na+nb)'],[1,theata2(1:na)'],Ts,'Variable','z^-1');
end
