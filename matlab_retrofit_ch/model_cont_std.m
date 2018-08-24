classdef model_cont_std < model
    properties
        params % Identify Number of Parameter 
        G      % Identification Model
        dG
        theta;
        n % denominator dimension (model dimension)
        m % numerator-denominator dimension
        nth
        has_x0
        x0
        lsim_type = 'zoh'; % Lsim type
    end
    
    methods
        function obj= model_cont_std(n, m, has_x0)
            if nargin < 3
                has_x0 = false;
            end
            obj.has_x0 = has_x0;
            a = sym('a',[n,1]);
            b = [sym('b0');sym('b',[n-m,1])];
            p = sym('p');
            obj.params = cellfun(@char,num2cell([a;b]),'UniformOutput',false);
            obj.nth = numel(obj.params);
            den = (p.^(n:-1:0))*[1;a];
            num = (p.^(n-m:-1:0))*b;
            
            obj.G = num/den;
            obj.theta = zeros(obj.nth, 1);
            pp = 0.1;
            eq = p+pp;
            for itr = 1:n-1
                pp = 1;
                %                obj.theta(itr) = nchoosek(n, itr);
                eq = eq*(p+pp);
            end
            eq = expand(eq);
            d = coeffs(eq);
            for itr = 1:n
                obj.theta(itr) = d(end-itr);
            end
            obj.theta(end) = 1e-5;
            obj.n = n;
            obj.m = m;
            for itr = 1:obj.nth
                obj.dG{itr} =  diff(obj.G,obj.params{itr});
            end
            if has_x0
                params_x0 = cell(n, 1);
                for itr = 1:n
                    params_x0{itr} = strcat('x0_', num2str(itr));
                end
                obj.params = [obj.params; params_x0];
            end
            obj.x0 = zeros(n, 1);
            
        end
        
        function [yhat, dyhat] = sim(obj, t, u, theta)
            %             num=numel(obj.params);
            % m=size(us,2);
            % dG = tf(zeros(n,m));
            den = [1, theta(1:obj.n)'];
            num = theta(obj.n+1:obj.nth)';
            sys = tf(num, den);
            yhat = lsim(sys,u,t, obj.lsim_type);
            dyhat = zeros(numel(t),numel(obj.params));
            
            for k=obj.n+1:obj.nth
                num = num*0;
                num(k-obj.n) = 1;
                dsys = tf(num, den);
                dyhat(:,k) =  lsim(dsys,u,t, obj.lsim_type);
            end
            for k=1:obj.n
                num = den*0;
                num(k+1) = 1;
                sys2 = tf(num, den);
                dsys = -sys*sys2;
                dyhat(:,k) = lsim(dsys,u,t, obj.lsim_type);
            end
            if obj.has_x0
                num = theta(obj.nth+1:end)';
                sys_init = tf(num, den);
                yhat_init = impulse(sys_init, t);
                yhat = yhat + yhat_init;
                
                for k = 1:obj.n
                    num =  num*0;
                    num(k) = 1;
                    dsys_init = tf(num, den);
                    dyhat(:, obj.nth+k) = impulse(dsys_init, t);
                end
                
                for k=1:obj.n
                    num = den*0;
                    num(k+1) = 1;
                    sys2 = tf(num, den);
                    dsys = -sys_init*sys2;
                    dyhat(:,k) = dyhat(:, k) + impulse(dsys,t);
                end
                
            end
            
        end
        
        function theta = get_params(obj)
            theta = obj.theta;
            if obj.has_x0
                theta = [theta; obj.x0];
            end
        end
        
        function set_params(obj, theta)
            nth = numel(obj.theta);
            obj.theta = theta(1:nth);
            if obj.has_x0
                obj.x0 = theta(nth+1:end);
            end
        end
        
        function Gtf = get_tf(obj, theta)
            if nargin == 1 || isempty(theta)
                theta = obj.get_params();
            else
                theta = theta(:);
            end
            %             Gtf = tools.sym2tf(obj.G, obj.params, obj.get_params());
            den = [1, theta(1:obj.n)'];
            num = theta(obj.n+1:obj.nth)';
            Gtf = tf(num, den);
        end
        
        function sys = get_sys(obj, varargin)
            sys = obj.get_tf(varargin{:});
        end
        
        function dG = get_tf_diff(obj)
            dG = tf(ones(numel(obj.dG) ,1));
            for itr = 1:numel(obj.dG)
                dG(itr) = tools.sym2tf(obj.dG{itr}, obj.params, obj.get_params());
            end
        end
        
        function fit(obj, t, u, y, theta)
            if nargin < 5
                if isempty(obj.theta) || any(isnan(obj.theta))
                    theta = obj.get_init(t, u, y);
                else
                    theta = obj.get_params();
                end
            end
            obj.fit@model(t, u, y, theta);
        end
        
        function fit_ls(obj, t, u, y, theta)
            if nargin < 5
                if isempty(obj.theta) || any(isnan(obj.theta))
                    theta = obj.get_init(t, u, y);
                else
                    theta = obj.get_params();
                end
            end
            obj.fit_ls@model(t, u, y, theta);
        end
        
        function theta = get_init(obj, t, u, y)
            [num,den]= tfdata(arx(iddata(y,u,mode(diff(t))),[obj.n,obj.n-obj.m+1,0]),'v');
            G_init = tf(num,den,mode(diff(t)));
            [num,den] = tfdata(minreal(d2c(G_init,'tustin')),'v');
            a_init = den(2:end);
            b_init = num(obj.m+1:end);
            theta = [a_init(:);zeros(obj.n-numel(a_init),1);b_init(:);zeros(obj.n-obj.m+1-numel(b_init),1);];
            
        end
        
        function set_tf(obj, G)
            if obj.n ~= order(G)
                %                 error('Model order is invalid.')
                n = obj.n;
                m = obj.m;
                [num, den] = tfdata(G, 'v');
                num = [num, 0.01*ones(1, 1+n-numel(num))];
                den = [den, zeros(1, 1+n-numel(den))];
                obj.theta = [den(2:end)/den(1), num(obj.m+1:end)/den(1)]';
                
            else
                [num, den] = tfdata(G, 'v');
                obj.theta = [den(2:end)/den(1), num(obj.m+1:end)/den(1)]';
            end
        end
    end
    
end

