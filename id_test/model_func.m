function cost = model_func(params,dim,u,y,t,Ts,lsim_type)
%{
load('iddata')
load('initial_tf')
params = [den(2:end),num];
u = v0;
y = w0;
N = 50000;
Ts = 0.25;
t = (0:N-1)'*Ts;
dim=6;
lsim_type='foh';
%}
    den = [1,params(1:dim)];
    num = params(dim+1:end);
    
    sys = tf(num,den);
    
    yhat = lsim(sys,u,t,lsim_type);
    %fpass = 0.1/(2*pi);
    %yhat = lowpass(yhat,fpass,Ts);
    
    error = y-yhat;
    %cost = error;
    cost = norm(error,2)/length(u);
    %cost = sum(error);
    
    %{
    dyhat = zeros(length(t),numel(params));
    for k=dim+1:numel(params)
        num = num*0;
        num(k-dim) = 1;
        dsys = tf(num, den);
        dyhat(:,k) =  lsim(dsys,u,t,lsim_type);
    end
    for k=1:dim
        num = den*0;
        num(k+1) = 1;
        sys2 = tf(num, den);
        dsys = -sys*sys2;
        dyhat(:,k) = lsim(dsys,u,t,lsim_type);
    end
    
    cost = 0;
    error = [yhat,dyhat]-y;
    for k = 1 : numel(params)+1
        cost = cost + norm(error(:,k),2);
    end
    %}
end