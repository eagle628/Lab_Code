classdef model < handle
    %MODEL
    
    properties(Abstract)
        params
    end
    
    properties
        fixed_params = {};
        fixed_values = [];
        str_display = 'iter';
        method = 'trust-region'; % obj.method;
        check_dv = 'off';
        gamma = 0;
        max_iter = 50000;
        np = 1000;
        lambda_mqt = 1;
        outfun=[];
        out_mqt;
        %         sim_p = @(obj, t, u) obj.sim(t, u, obj.get_params());
    end
    
    methods(Abstract)
        sim(obj)
        set_params(obj, params)
        get_params(obj)
    end
    methods
        
        function varargout = sim_p(obj, t, u)
            varargout = cell(nargout, 1);
            [varargout{:}] = obj.sim(t, u, obj.get_params());
        end
        
        function varargout = sim_fix_p(obj, t, u)
            varargout = cell(nargout, 1);
            [varargout{:}] = obj.sim_fix(t, u, obj.get_params_fixed());
        end
        
        function [val, dval] = eval_func(obj, t, u, y, theta, weight)
            if nargin < 6
                weight = 1;
            end
            try
            [yhat, dyhat] = obj.sim_fix(t, u, theta);
            e = (y(:) - yhat(:));
            val = e' * (weight(:).* e);
            dval = -2*dyhat'*(weight(:).*e);
            catch
                val = inf;
                dval = nan;
            end
            
        end
        
        function [theta, fval] = fit(obj, t, u, y, theta0)
            if nargin < 5
                theta0 = obj.get_params_fixed();
            else
                theta0 = obj.get_params_fixed(theta0);
            end
            func = @(x) obj.eval_func(t, u, y, x);
            
            options = optimoptions(@fminunc,'GradObj','on',...
                'Display', obj.str_display, 'Algorithm', obj.method,...
                'TolFun', 1e-8, 'tolX', 1e-8,...
                'MaxFunEvals', 1e4, 'MaxIter', 1e4,...
                'DerivativeCheck', obj.check_dv,...
                'FinDiffRelStep', 1e-5,...
                'OutputFcn', @(x, opt, state) callback_fmin(x, opt, state, y(:), @(x) obj.sim_fix(t, u, x), obj)...
                );
            [theta, fval] = fminunc(func, theta0, options);
            obj.set_params_fixed(theta);
        end
        
        function [theta, fval] = fit_lbfgs(obj, t, u, y, theta0, lb, ub)
            if nargin < 5 || isempty(theta0)
                theta0 = obj.get_params_fixed();
            else
                theta0 = obj.get_params_fixed(theta0);
            end
            if nargin<6
                lb = -inf(size(theta0));
                ub = inf(size(theta0));
            end
            func = @(x) obj.eval_func(t, u, y, x);
            opts.x0 = theta0;
            opts.maxIts=1e6;
            opts.maxTotalIts=1e6;
            opts.pgtol = 1e-20;
            opts.factr = 1e-15;
            opts.printEvery = 200;
            [theta, fval, info] = lbfgsb(func, lb, ub, opts);
            %             options = optimoptions(@fminunc,'GradObj','on',...
            %                 'Display', obj.str_display, 'Algorithm', obj.method,...
            %                 'TolFun', 1e-30, 'tolX', 1e-20,...
            %                 'MaxFunEvals', 1e4, 'MaxIter', 1e4,...
            %                 'DerivativeCheck', obj.check_dv,...
            %                 'FinDiffRelStep', 1e-5);
            %             [theta, fval] = fminunc(func, theta0, options);
            obj.set_params_fixed(theta);
        end
        
        function [theta, fval] = fit_con(obj, t, u, y, theta0, lb, ub)
            if nargin < 5 || isempty(theta0)
                theta0 = obj.get_params_fixed();
            else
                theta0 = obj.get_params_fixed(theta0);
            end
            if nargin<6
                lb = -inf(size(theta0));
                ub = inf(size(theta0));
            end
            func = @(x) obj.eval_func(t, u, y, x);
            options = optimoptions(@fmincon,'GradObj','on',...
                'Display', obj.str_display, ...
                'TolFun', 1e-30, 'tolX', 1e-20,...
                'MaxFunEvals', 1e4, 'MaxIter', 1e4,...
                'DerivativeCheck', obj.check_dv,...
                'FinDiffRelStep', 1e-5,...
                'OutputFcn', @(x, opt, state) callback_fmin(x, opt, state, y(:), @(x) obj.sim(t, u, x), obj)...
                );
            [theta, fval] = fmincon(func, theta0, [], [], [], [],lb, ub,[], options);
            obj.set_params_fixed(theta);
            
        end
        
        function theta = fit_ls(obj, t, u, y, theta0)
            if nargin < 5
                theta0 = obj.get_params_fixed();
            else
                theta0 = obj.get_params_fixed(theta0);
            end
            func = @(x) obj.sim_ls(t, u, y, x);
            options = optimoptions(@lsqnonlin, 'Jacobian','on',...
                'Display', obj.str_display, 'Algorithm', 'levenberg-marquardt',...
                'DerivativeCheck', obj.check_dv,...
                'FinDiffRelStep', 1e-5,...
                'OutputFcn', @(x, opt, state) callback_fmin(x, opt, state, y(:), @(x) obj.sim(t, u, x), obj)...
                );
            theta = lsqnonlin(func, theta0, [], [], options);
            obj.set_params_fixed(theta);
        end
        
        [theta ,Jhistory] = fit_sgd(obj, t, u, y, theta, learning_ratio, weight);
        [theta ,Jhistory] = fit_momentum(obj, t, u, y, theta, learning_ratio, mu, weight);
        [theta ,Jhistory] = fit_adam(obj, t, u, y, theta, learning_ratio, rho1, rho2, epsilon, weight, snr);
        [theta ,Jhistory] = fit_adamax(obj, t, u, y, theta, learning_ratio, rho1, rho2, weight, snr);
        [theta, Jhistory] = fit_Santa_E(obj, t, u, y, theta, learning_ratio, rho, epsilon, burnin, beta_func, weight_b);
        [theta, Jhistory] = fit_Santa_S(obj, t, u, y, theta, learning_ratio, rho, epsilon, burnin, beta_func, weight_b);
        
        function theta = fit_mymqt(obj, t, u, y, theta0, weight)
            %             if nargin < 5
            %                 if any(isnan(obj.get_params()))
            %                     theta0 = obj.params_default;
            %                 else
            %                     theta0 = obj.get_params();
            %                 end
            %             end
            
            if nargin < 5 || isempty(theta0)
                [theta0,p] = obj.get_params_fixed();
            else
                [theta0,p] = obj.get_params_fixed(theta0);
            end
            model = struct();
            if nargin == 6
                model.w = weight;
            end
            
            model.params = p;
            model.params_init = theta0;
            model.y = y(:);
            model.sim = @(var, val) obj.sim_fix(t, u, val);
            options = struct();
            options.lambda = obj.lambda_mqt;
            options.maxitr = obj.max_iter;
            if numel(obj.gamma) == 1
                obj.gamma = ones(numel(obj.params), 1)*obj.gamma;
            end
            options.gamma = obj.get_params_fixed(obj.gamma);
            if strcmp(obj.str_display, 'none')
                options.callback = @(varargin) [];
            elseif ~isempty(obj.outfun)
                options.callback = obj.outfun;
            end
            out = est_param_mqt(model, options);
            theta = out.theta;
            obj.set_params_fixed(theta);
            obj.out_mqt = out;
        end
        
       
        function [e, de] = sim_ls(obj, t, u, y, theta)
            [yhat, dyhat] = obj.sim_fix(t, u, theta);
            if std(y) ~= 0
                d = std(y);
            else
                d = 1;
            end
            e = (y-yhat)/d;
            de = -dyhat/d;
        end
        
        function [y, dy] = sim_fix(obj, t, u, theta)
            theta_in = zeros(numel(obj.params), 1);
            [~, idx] = setdiff(obj.params, obj.fixed_params);
            idx = sort(idx);
            theta_in(idx) = theta;
            for itr = 1:numel(obj.fixed_params)
                theta_in(strcmp(obj.params, obj.fixed_params{itr})) = obj.fixed_values(itr);
            end
            if nargout == 1
                y = obj.sim(t, u, theta_in);
            else
            [ytmp, dytmp] = obj.sim(t, u, theta_in);
            l = size(ytmp, 2);
            n = numel(theta_in);
            N = size(ytmp, 1);
            dy = zeros(numel(ytmp), n);
            for itr = 1:l
               dy((1:N)+N*(itr-1),:) = dytmp(:, itr:l:end); 
            end            
            dy = dy(:, idx);
            y = ytmp(:);
            end
        end
        
        function add_fixed_params(obj, params, vals)
            if isempty(params)
                return;
            end
            %             obj.fixed_params = [obj.fixed_params; params];
            %             obj.fixed_values = [obj.fixed_values; vals];
            theta = obj.get_params();
            if iscell(params)
                for itr = 1:numel(params)
                    p = params{itr};
                    if ismember(p, obj.fixed_params)
                        obj.fixed_values(strcmp(p, obj.fixed_params)) = vals(itr);
                    else
                        obj.fixed_params = [obj.fixed_params; {p}];
                        obj.fixed_values = [obj.fixed_values; vals(itr)];
                    end
                    theta(strcmp(p, obj.params)) = vals(itr);
                end
            else
                if ismember(params, obj.fixed_params)
                    obj.fixed_values(strcmp(params, obj.fixed_params)) = vals;
                else
                    obj.fixed_params = [obj.fixed_params; {params}];
                    obj.fixed_values = [obj.fixed_values; vals];
                end
                theta(strcmp(params, obj.params)) = vals;
            end
            obj.set_params(theta);
        end
        
        function remove_fixed_params(obj, params)
            [~, idx] = intersect(obj.fixed_params, params);
            obj.fixed_params(idx) = [];
            obj.fixed_values(idx) = [];
        end
        
        function [theta, p] = get_params_fixed(obj, theta)
            if nargin < 2
                theta = obj.get_params();
            end
            [~, idx] = setdiff(obj.params, obj.fixed_params);
            idx = sort(idx);
            theta = theta(idx);
            p = obj.params(idx);
        end
        
        function theta_out = get_params_all(obj, theta)
            if nargin < 2
                theta_out = obj.get_params();
            else
                theta_out = zeros(numel(obj.params), 1);
                [~, idx] = setdiff(obj.params, obj.fixed_params);
                idx = sort(idx);
                theta_out(idx) = theta;
                for itr = 1:numel(obj.fixed_params)
                    theta_out(strcmp(obj.params, obj.fixed_params{itr})) = obj.fixed_values(itr);
                end
            end
        end
        
        function set_params_fixed(obj, theta)
            theta_in = obj.get_params_all(theta);
            obj.set_params(theta_in);
        end
        
        
        function theta = fit_pso(obj, t, u, y, theta0)
            if nargin < 5
                theta0 = obj.get_params_fixed();
            else
                theta0 = obj.get_params_fixed(theta0);
            end
            
            %             func = @(x) obj.fit_func(t, u, y, x);
            func = @(x) obj.eval_func(t, u, y, x);
            X = bsxfun(@plus,randn(obj.np, numel(theta0)),theta0');
            X(1,:) = theta0';
            theta = mypso_solve(func, X, obj.max_iter);
            obj.set_params_fixed(theta);
            
        end
        
        function out = get_sys(obj, varargin)
           out =  [];
        end
        
        function out = get_tf(obj, varargin)
           out = tf(obj.get_sys(varargin{:}));
        end
        
        function out = tf(obj)
            out = obj.get_tf();
        end
        
        function out = ss(obj)
           out = ss(obj.get_sys()); 
        end
        
        function set_tf(obj, sys)
           obj.set_sys(ss(sys)); 
        end
        
        function set_sys(obj, sys)
           []; 
        end
       
    end
    
    
end


