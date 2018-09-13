 classdef model_ss < model
    
    properties
        params
        gen_ss
        lsim_type = 'foh';
    end
    
    methods
        function obj = model_ss(gen_ss)
            obj.gen_ss = gen_ss;
            obj.params = gen_ss.params;
        end
        
        function set_sys(obj, sys)
            obj.gen_ss.set_sys(sys);
        end
        
        function n4sid(obj, t, u, y)
            z = iddata(y, u, mode(diff(t)));
            m = n4sid(z, order(ss(obj)));
            sys = d2c(m);
            if order(sys) > order(ss(obj))
                sys = balred(sys, order(ss(obj)));
            end
            obj.set_sys(sys);
        end
        
        function theta = get_params(obj)
            theta = obj.gen_ss.get_params();
        end
        
        function set_params(obj, theta)
            obj.gen_ss.set_params(theta);
        end
        
        function sys = get_sys(obj, theta)
            if nargin == 1 || isempty(theta)
                theta = obj.get_params;
            end
            [A, B, C, D] = obj.gen_ss.get_ss(theta);
            sys = ss(A, B, C, D);
        end
        
        function [y, dy] = sim(obj, t, u, theta)
            np = numel(obj.params);
            if nargout == 1
                [A, B, C, D] = obj.gen_ss.get_ss(theta);
                y = lsim(ss(A, B, C, D), u, t, obj.lsim_type);
            else
                [A, B, C, D, dA, dB, dC, dD] = obj.gen_ss.get_ss(theta);
                n = size(A, 1);
                l = size(C, 1);
                Abig = kron(eye(np+1), A);
                Bbig = kron(ones(np+1, 1), B*0);
                Bbig(1:n,:) = B;
                Cbig = kron(eye(np+1), C);
                Dbig = kron(ones(np+1,1), D*0);
                Dbig(1:l, :) = D;
                for itr  = 1:numel(obj.params)
                    Abig(itr*n+(1:n), 1:n) = dA{itr};
                    Bbig(itr*n+(1:n), :) = dB{itr};
                    Cbig(itr*l+(1:l), 1:n) = dC{itr};
                    Dbig(itr*l+(1:l), :) = dD{itr};
                end
                %             sys = ss(Abig, Bbig, Cbig, Dbig);
                %                 sys = minreal(ss(Abig, Bbig, Cbig, Dbig), 1e-8, false);
                sys = ss(Abig, Bbig, Cbig, Dbig);
                %             sys = minreal(ss(Abig, Bbig, Cbig, Dbig), 1e-10, false);
                y = lsim(sys, u, t, obj.lsim_type);
                dy = y(:, l+1:end);
                y = y(:, 1:l);
            end
        end
        
        function sys = get_sys_big(obj, theta)
            if nargin == 1 || isempty(theta)
                theta = obj.get_params();
            end
            np = numel(theta);
            [A, B, C, D, dA, dB, dC, dD] = obj.gen_ss.get_ss(theta);
            n = size(A, 1);
            l = size(C, 1);
            Abig = kron(eye(np+1), A);
            Bbig = kron(ones(np+1, 1), B*0);
            Bbig(1:n,:) = B;
            Cbig = kron(eye(np+1), C);
            Dbig = kron(ones(np+1,1), D*0);
            Dbig(1:l, :) = D;
            for itr  = 1:numel(obj.params)
                Abig(itr*n+(1:n), 1:n) = dA{itr};
                Bbig(itr*n+(1:n), :) = dB{itr};
                Cbig(itr*l+(1:l), 1:n) = dC{itr};
                Dbig(itr*l+(1:l), :) = dD{itr};
            end
            sys = ss(Abig, Bbig, Cbig, Dbig);
            sys.OutputGroup.y = 1:size(C, 1);
            sys.OutputGroup.dy = size(C, 1)+1:size(Cbig, 1);
        end
        
        function fit_balance(obj, t, u, y, theta0)
            if nargin < 5
                theta0 = obj.get_params_fixed();
            else
                theta0 = obj.get_params_fixed(theta0);
            end
            func = @(x) obj.eval_func(t, u, y, x);
            
            options = optimoptions(@fmincon,'SpecifyConstraintGradient',true,...
                'SpecifyObjectiveGradient',true,...
                'Algorithm', 'interior-point', ...
                'HessianApproximation', 'bfgs', ...
                'Display', obj.str_display,...
                'OptimalityTolerance', 1e-12, ...
                'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4,...
                'DerivativeCheck', obj.check_dv,...
                'ConstraintTolerance', 1e-3,...
                'UseParallel', true,...
                'OutputFcn', @(x, opt, state) callback_fmin(x, opt, state, y(:), @(x) obj.sim(t, u, x), obj)...
                );
            [theta, fval] = fmincon(func, theta0, [], [], [], [], [], [], @(x) obj.gen_ss.con_blance(x), options);
            obj.set_params_fixed(theta);
        end
        
        function fit_constraint(obj, t, u, y, theta0)
            if nargin < 5
                theta0 = obj.get_params_fixed();
            else
                theta0 = obj.get_params_fixed(theta0);
            end
            func = @(x) obj.eval_func(t, u, y, x);
            
            options = optimoptions(@fmincon,'SpecifyConstraintGradient',true,...
                'SpecifyObjectiveGradient',true,...
                'Algorithm', 'interior-point', ...
                'HessianApproximation', 'bfgs', ...
                'Display', obj.str_display,...
                'OptimalityTolerance', 1e-12, ...
                'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4,...
                'DerivativeCheck', obj.check_dv,...
                'ConstraintTolerance', 1e-3,...
                'UseParallel', true,...
                'StepTolerance', 1e-12, ... %                 'SubproblemAlgorithm','cg',...
                'OutputFcn', @(x, opt, state) callback_fmin(x, opt, state, y(:), @(x) obj.sim(t, u, x), obj)...
                );
            [theta, fval] = fmincon(func, theta0, [], [], [], [], [], [], @(x) obj.gen_ss.constraint(x), options);
            obj.set_params_fixed(theta);
        end
        
        function find_stable(obj, t, u, y, theta0)
            if nargin < 5
                theta0 = obj.get_params_fixed();
            else
                theta0 = obj.get_params_fixed(theta0);
            end
            func = @(x) obj.eval_dummy(t, u, y, x);
            
            options = optimoptions(@fmincon,'SpecifyConstraintGradient', false,...
                'SpecifyObjectiveGradient',true,...
                'Algorithm', 'interior-point', ...
                'HessianApproximation', 'bfgs', ...
                'Display', obj.str_display,...
                'OptimalityTolerance', 1e-12, ...
                'MaxFunctionEvaluations', 1e4, 'MaxIterations', 1e4,...
                'DerivativeCheck', obj.check_dv,...
                'ConstraintTolerance', 1e-6,...
                'UseParallel', true,...
                'OutputFcn', @(x, opt, state) callback_fmin(x, opt, state, y(:), @(x) obj.sim(t, u, x), obj)...
                );
            [theta, fval] = fmincon(func, theta0, [], [], [], [], [], [], @(x) obj.gen_ss.con_stable(x), options);
            obj.set_params_fixed(theta);
        end
        
        function [v, dv] = eval_dummy(obj, varargin)
            v = 0;
            dv = zeros(numel(obj.get_params_fixed), 1);
        end
        
        
        
        function [val, dval] = eval_func(obj, t, u, y, theta, weight)
            if nargin < 6
                weight = 1;
            end
            try
                %             try
                %                 [yhat, dyhat] = obj.sim_fix(t, u, theta);
                %                 e = (y(:) - yhat(:));
                theta = obj.get_params_all(theta);
                
                [A, B, C, D, dA, dB, dC, dD] = obj.gen_ss.get_ss(theta); % So rectifier obj's get_ss, A dim = local+env
                Ts = mode(diff(t));
                nx = size(A, 1);
                nu = size(B, 2);
                ny = size(C, 1);
                [~, idx] = setdiff(obj.params, obj.fixed_params);
                idx = sort(idx(:));
                
                
                switch(obj.lsim_type)
                    case 'foh'
                        M1 =  [A , B , zeros(nx,nu)  ; ...
                            sparse(nu,nx+nu)  speye(nu)/Ts ; ...
                            sparse(nu,nx+2*nu)];
                        Md = expm(M1*Ts);
                        F1 = Md(1:nx,nx+1:nx+nu);
                        F2 = Md(1:nx,nx+nu+1:nx+2*nu);
                        
                        Ad = Md(1:nx,1:nx);
                        Bd = F1 + Ad*F2 - F2;
                        Cd = C;
                        Dd = D + C*F2;
                        yxhat = lsim(ss(Ad, Bd, [Cd; eye(size(A))], [Dd; tools.zeros(A, D)], Ts), u);
                        yhat = yxhat(:, 1:ny);
                        xhat = yxhat(:, ny+1:end);
                        
                        e = (y - yhat);
                        val = sum(sum(weight.*(e.^2)));
                        
                        if nargout < 2
                            return
                        end
                        sys_backward = ss(Ad', Cd', eye(size(Ad)), 0, Ts);
                        Xbar = flipud(lsim(sys_backward, flipud(weight.*e)));
                        dval = zeros(numel(idx), 1);
                        for itr = 1:numel(idx)
                            i = idx(itr);
                            if all(all(dA{i}==0)) && all(all(dB{i}==0))
                                dAi = sparse(nx, nx);
                                dBi = sparse(nx, nu);
                                dDi = sparse(dD{i} + dC{i}*F2);
                            else
                                M1 =  [A , B , zeros(nx,nu)  ; ...
                                    sparse(nu,nx+nu)  speye(nu)/Ts ; ...
                                    sparse(nu,nx+2*nu)];
                                M2 =  [dA{i} , dB{i} , sparse(nx,nu)  ; ...
                                    sparse(nu,nx+nu)  sparse(nu, nu); ...
                                    sparse(nu,nx+2*nu)];
                                M = [M1, M2; tools.zeros(M2, M1), M1];
                                Md = expm(M*Ts);
                                F1 = Md(1:nx,nx+1:nx+nu);
                                F2 = Md(1:nx,nx+nu+1:nx+2*nu);
                                
                                Ad = Md(1:nx,1:nx);
                                %                     Bd = F1 + Ad*F2 - F2;
                                %                     Cd = C;
                                %                     Dd = D + C*F2;
                                dAi = sparse(Md(1:nx, (nx+2*nu)+(1:nx)));
                                dF1 = sparse(Md(1:nx, 2*nx+2*nu+(1:nu)));
                                dF2 = sparse(Md(1:nx, 2*nx+3*nu+(1:nu)));
                                dBi = sparse(dF1-dF2+Ad*dF2+dAi*F2);
                                dDi = sparse(dD{i}+dC{i}*F2+C*dF2);
                            end
                            dCi = sparse(dC{i});
                            W = xhat*dAi' + u*dBi';
                            a = sum(sum(W.*Xbar));
                            b = sum(sum((xhat*dCi'+u*dDi').*weight.*e));
                            dval(itr) = -2*(a+b);
                        end
                        %                     M =  [A , B , zeros(nx,nu)  ; ...
                        %                         zeros(nu,nx+nu)  eye(nu)/Ts ; ...
                        %                         zeros(nu,nx+2*nu)];
                        %                     Md = expm(M*Ts);
                        %                     F1 = Md(1:nx,nx+1:nx+nu);
                        %                     F2 = Md(1:nx,nx+nu+1:nx+2*nu);
                        %
                        %                     Ad = Md(1:nx,1:nx);
                        %                     Bd = F1 + Ad*F2 - F2;
                        %                     Cd = C;
                        %                     Dd = D + C*F2;
                        
                    case 'zoh'
                        Md = expm([A, B;sparse(nu, nx+nu)]*Ts);
                        Ad = Md(1:nx, 1:nx);
                        Bd = Md(1:nx, nx+(1:nu));
                        Cd = C;
                        Dd = D;
                        yxhat = lsim(ss(Ad, Bd, [Cd; eye(size(A))], [Dd; tools.zeros(A, D)], Ts), u);
                        yhat = yxhat(:, 1:ny);
                        xhat = yxhat(:, ny+1:end);
                        
                        e = (y - yhat);
                        val = sum(sum(weight.*(e.^2)));
                        
                        if nargout < 2
                            return
                        end
                        sys_backward = ss(Ad', Cd', eye(size(Ad)), 0, Ts);
                        Xbar = flipud(lsim(sys_backward, flipud(weight.*e)));
                        dval = zeros(numel(theta), 1);
                        for itr = 1:numel(idx)
                            i = idx(itr);
                            if all(all(dA{i}==0)) && all(all(dB{i}==0))
                                dAi = sparse(nx, nx);
                                dBi = sparse(nx, nu);
                            else
                                Mbig = [A, B, dA{i}, dB{i};sparse(nu, 2*(nx+nu));
                                    sparse(nx, (nx+nu)), A, B; sparse(nu, 2*(nx+nu))];
                                Md = expm(Mbig*Ts);
                                dAi = sparse(Md(1:nx, (nx+nu)+(1:nx)));
                                dBi = sparse(Md(1:nx, (nx+nu)+nx+(1:nu)));
                            end
                            dCi = sparse(dC{i});
                            dDi = sparse(dD{i});
                            W = xhat*dAi' + u*dBi';
                            a = sum(sum(W.*Xbar));
                            b = sum(sum((xhat*dCi'+u*dDi').*weight.*e));
                            dval(itr) = -2*(a+b);
                        end
                end
            catch
                val = inf;
                dval = nan;
            end
        end
        
    end
    
end



