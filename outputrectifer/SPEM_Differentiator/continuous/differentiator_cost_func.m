function [cost,gra] = differentiator_cost_func(dim,params,data,K,Filter,sym_gra)
warning('off')
    if isempty(Filter)
        Filter = tf(1);
    end
    den = [1,params(1:dim)];
    num = params(dim+1:end);
    num = [num,0];
    
    P_hat = tf(num,den);

    S_hat = 1/(1+P_hat*K);
    G_eu = -Filter*S_hat*P_hat;
    G_ey = Filter*S_hat;

    G_eu = minreal(G_eu, [], false);
    G_ey = minreal(G_ey, [], false);

    time = 0:data.Ts:data.Ts*(data.N-1);
    e = lsim(G_eu,data.u,time)+lsim(G_ey,data.y,time);

    cost = e./data.N;
    if nargout == 2
        nvar = size(sym_gra,1);
        gra = zeros(size(e));
        for itr = 1 : nvar
            dG_ey = Filter*sym2tf(sym_gra{itr,1}, params);
            dG_ey = minreal(dG_ey, [], false);
            dG_eu = Filter*sym2tf(sym_gra{itr,2}, params);
            dG_eu = minreal(dG_eu, [], false);
            gra(:,itr) = lsim(dG_eu,data.u,time) + lsim(dG_ey,data.y,time);
        end
%             gra = mean(gra,1);
    end
warning('on')
end
