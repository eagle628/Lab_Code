function [cost, gra] = SPEM_cost_func(params,data,sys_rect,dim)
    den = [1,params(1:dim)];
    num = params(dim+1:end);
    B = get_tf(num,data.Ts);
    F = get_tf(den,data.Ts);
    
    Ge = B/F;
    %Ge = tf(num,den,data.Ts);
    Gl = sys_rect({'w'});
    %Gl = 0.1;
    Filter = sys_rect({'omega'});
%     Filter = tf(1);
    
    G1 =   1 / ( 1 - Ge*Gl);
    G2 = -Ge / ( 1 - Gl*Ge);
    
    out_x1 = lsim(Filter*G1,data.y,'foh');
    out_x2 = lsim(Filter*G2,data.u,'foh');
    out_x = out_x1 + out_x2;
    
%     cost = norm(out_x)/data.N;
    cost = out_x/data.N;

    % Gradient Or Jacbian
%     if nargout > 1
%         gra = zeros(numel(params),1);
%         idx = 1;
%         for itr = 1 : dim
%             dF_u = -B/(F*F)*get_tf([zeros(1,itr),1],data.Ts);
%             gra(idx) = mean(lsim( dF_u, data.u));
%             idx = idx + 1;
%         end
%         for itr = 1 : numel(params)-dim
%             dB_u = 1/F*get_tf([zeros(1,itr-1),1],data.Ts);
%             gra(idx) = mean(lsim( dB_u, data.u));
%             idx = idx + 1;
%         end
%     end
end
%% Local
function G = get_tf(param,Ts)
    G = tf(param,1,Ts,'Variable','q^-1');
end
