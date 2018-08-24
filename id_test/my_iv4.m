%function [A,B] = my_iv4(u,y,system_iv,Ts)
clear all
load('iddata');
y = v0;
u = w0;
system_iv = [6,7,0];
Ts = 0.25;

    % preparation
    number = length(u);
    t = 0:Ts:Ts*(number-1);
    na = system_iv(1);
    nb = system_iv(2);
    nk = system_iv(3);
    %% Stage I
    %system_iv = [na,nb,nk];
    G1 = myARX(u,y,[na,nb,nk],Ts);
    
    %% Stage II
    x1 = lsim(G1,u,t);
    d = nk-1;
    X = zeros(number,na);
    Y = zeros(number,na);
    U = zeros(number,nb);
    for i = 1 : na
        X(:,i) = [zeros(i,1);y(1:end-i)];
    end
    for i = 1 : na
        Y(:,i) = [zeros(i,1);y(1:end-i)];
    end
    for i = 1 : nb
        U(:,i) = [zeros(i+d,1);u(1:end-d-i)];
    end
    zeta1 = [-X,U];
    PHI = [-Y,U];
    
    theta2 = (zeta1'*PHI)\(zeta1'*y);
    G2 = tf([zeros(1,nk),theta2(na+1:na+nb)'],[1,theta2(1:na)'],Ts,'Variable','z^-1');
    %% Step III
    A2 = G2.Denominator{:};
    B2 = G2.Numerator{:};
    
    Y_prime = [y,Y];
    %U_prime = [u,U];
    
    w2 = A2*Y_prime' - B2*U';
    w2 = w2';
    G_high = myARX(u,y,[na+nb,1,nk],Ts);
    x2 = lsim(G2,u,t);
    x2_h = lsim(G_high,u,t);
    e = x2-x2_h;
    
    W = zeros(number,na+nb);
    for i = 1 : na+nb
        W(:,i) = [zeros(i,1);w2(1:end-i)];
    end
    L = W\e;
    %% Step IV
    X2 = zeros(number,na);
    for i = 1 : na
        X2(:,i) = [zeros(i,1);x2(1:end-i)];
    end
    zeta2 = L'*[-X2,U]';
    zeta2 = zeta2';
    
    PHI_f = L'*PHI';
%end

%% local function
function iv7_118()
end