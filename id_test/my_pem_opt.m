function [model] = my_pem_opt(u,y,system,Ts)
    t = (0:length(u)-1)'*Ts;
    na = system(1);
    %% initial condition By my ARX
    [D_sys] = myARX(u,y,system,Ts);
    sys = d2c(D_sys);
    %data = iddata(y,u,Ts);
    %sys = d2c(tf(oe(data,system))); 
    % これを初期値として，optimization toolboxでPEMを行う
    num = cell2mat(sys.Numerator);
    den = cell2mat(sys.Denominator);

    % optimiztion option
    opt = optimoptions('fminunc');
    opt.Display = 'iter-detailed';
    opt.FiniteDifferenceType = 'forward';
    opt.FunctionTolerance = 1e-6;
    opt.StepTolerance = 1e-6;
    opt.MaxFunctionEvaluations = 500*(numel(num)+numel(den));

    % Cost Functino Assign
    params_ini = [den(2:end),num];
    cost_func = @(params)model_func(params,numel(den)-1,u,y,t,Ts,'foh');
    
    [params,fval,exitflag,output] = fminunc(cost_func,params_ini,opt)

    den_opt = [1,params(1:na)];
    num_opt = params(na+1:end);

    model = tf(num_opt,den_opt);
    %model= minreal(ss(tf(num_opt,den_opt)));
end