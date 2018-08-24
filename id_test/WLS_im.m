clear all
close all


N = 50000;
Ts = 0.25;
t = (0:N-1)'*Ts;

load('iddata')
load('original')
load('ori_low')

%% pem
dim = 6;
na = dim;nb = dim+1;nc = dim;nd = dim;nf = dim;nk = 0;
nn = [na,nb,nc,nd,nf,nk];
data0 = iddata(v0,w0,Ts);
%model1 = pem(data0,nn,'Display','on');

%% initial condition By my ARX
system = [dim,dim+1,0];
[D_sys] = myARX(w0,v0,system,Ts);
sys = d2c(D_sys);
% これを初期値として，optimization toolboxでPEMを行う
%sys = tf(lower);
sys = d2c(tf(oe(data0,[6,6,0]))); 
num = cell2mat(sys.Numerator);
den = cell2mat(sys.Denominator);

% optimiztion option
opt = optimoptions('fminunc');
opt.Display = 'iter-detailed';
opt.FiniteDifferenceType = 'forward';
opt.FunctionTolerance = 1e-6;
opt.StepTolerance = 1e-6;
opt.MaxFunctionEvaluations = 500*(numel(num)+numel(den));

% Assign
params_ini = [den(2:end),num];
cost_func = @(params)model_func(params,dim,w0,v0,t,Ts,'foh');
params = fminunc(cost_func,params_ini,opt);

den_opt = [1,params(1:dim)];
num_opt = params(dim+1:end);

model2 = minreal(ss(tf(num_opt,den_opt)));

%% plot area
figure
bode(sys_env)
hold on
%bode(model1)
bode(model2)
bode(lower)

%% fault
%{
params = 13;
max_itr = 1000;
tolerance = 1e-3;
for itr = 1: max_itr
    [vhat, dvhat] = sim(t, w0, num,den,params);
    cost = norm(v0-vhat,2);
    Ja = v0-dvhat;
    dx = -(Ja'*Ja)\(Ja'*(v0-vhat));
    dx = dx';
    if norm(dx,inf) < tolerance
        break;
    end
    den = den + [0,dx(1:na)];
    num = num + dx(na+1:end);
        
end
%}
%% 
function [yhat, dyhat] = sim(t, u, num,den,params)
    model_n = numel(den)-1;
    lsim_type = 'foh';
    sys = tf(num, den);

    yhat = lsim(sys,u,t, lsim_type);
    dyhat = zeros(numel(t),numel(params));

    for k=model_n+1:params
        num = num*0;
        num(k-model_n) = 1;
        dsys = tf(num, den);
        dyhat(:,k) =  lsim(dsys,u,t,lsim_type);
    end
    for k=1:model_n
        num = den*0;
        num(k+1) = 1;
        sys2 = tf(num, den);
        dsys = -sys*sys2;
        dyhat(:,k) = lsim(dsys,u,t,lsim_type);
    end
end
%% function
function [G_func,params] = tf_jacobian(n,m)
    s = sym('s');
    a = sym('a',[n,1]);
    b = [sym('b0');sym('b',[n-m,1])];
    c = [a;b];
    
    num = (s.^(0:1:n-m))*b;    den = (s.^(0:1:n))*[1;a];

    G_sym = num/den;

    dG_sym = jacobian(G_sym,[a;b]);
    dG_sym = expand(dG_sym);
    dG_sym = collect(dG_sym,'s');
    params = length([a;b]);
    G_func = cell(params,2);
    for i = 1:params
         [num_sym,den_sym] = numden(dG_sym(i));
          num_coe = coeffs(num_sym,'s');
          den_coe = coeffs(den_sym,'s');
          G_func(i,1) = {symfun(num_coe,c)};
          G_func(i,2) = {symfun(den_coe,c)};
    end
end
