function cost = SPEM_cost_func2(params,data,sys_rect,dim)
    den = [1,params(1:dim)];
    num = params(dim+1:end);
    
    Ge = tf(num,den);
    Gl = sys_rect({'w'});
    Filter = sys_rect({'omega'});
%     Filter = tf(1);

    G1 =    1 / ( 1 - Ge*Gl);
    G2 =  -Ge / ( 1 - Gl*Ge);
    
    time = 0:data.Ts:data.Ts*(data.N-1);
    out_x1 = lsim( Filter*G1, data.y, time, 'foh');
    out_x2 = lsim( Filter*G2, data.u, time, 'foh');
    out_x = out_x1 + out_x2;
    
%     cost = norm(out_x)/data.N;
    cost = out_x/data.N;
end