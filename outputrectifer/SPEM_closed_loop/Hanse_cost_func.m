function [cost] = Hanse_cost_func(params,data,K,dim)
try
    den = [1,params(1:dim)];
    num = params(dim+1:end);
%     num(1) = 0;
    
    P_hat = tf(num,den,data.Ts);
    
    S_hat = 1/(1+P_hat*K);
    G_eu = -S_hat*P_hat;
    G_ey = S_hat;
    
    G_eu = minreal(G_eu, [], false);
    G_ey = minreal(G_ey, [], false);
    
    e = lsim(G_eu,data.u)+lsim(G_ey,data.y);
    
    
%     cost = norm(e)/data.N;
    cost = e./data.N;
    
catch
   cost = inf; 
end
    
    
end