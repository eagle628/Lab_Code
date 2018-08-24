function cost = out_rect_func(params,w,v,G_wv,G_yv,Ts,dim)
    den = [1,params(1:dim)];
    num = params(dim+1:end);
    
    Ge = tf(num,den,Ts);
    
    Gxv =  G_yv/(1-G_wv*Ge);
    Gxw = -G_yv*Ge/(1-G_wv*Ge);
    
    out_x = lsim(Gxv,v)+lsim(Gxw,w);
    
    cost = norm(out_x,2)/length(out_x);
    
end