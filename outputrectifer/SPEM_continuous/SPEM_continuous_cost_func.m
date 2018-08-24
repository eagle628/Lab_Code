function [cost,gra] = SPEM_continuous_cost_func(dim,params,data,K,Filter,sym_gra)
warning('off')
    if isempty(Filter)
        Filter = tf(1);
    end
    den = [1,params(1:dim)];
    num = params(dim+1:end);

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

%% Local
function tf_sys = sym2tf(sym_sys,params)
    var = symvar(sym_sys);
%     a = sym('a',[dim, 1]);
%     b = sym('b',[dim+1, 1]);
%     sym_sys = subs(sym_sys, [a ; b], params');
    sym_sys = subs(sym_sys, var(1:end-1), params');
    [num_s, den_s] = numden(sym_sys);
%     syms s
    num = double(fliplr(coeffs(num_s, var(end))));
    den = double(fliplr(coeffs(den_s, var(end))));
    tf_sys = minreal(tf(num, den), [], false);
end