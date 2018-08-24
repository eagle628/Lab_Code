classdef pem_discrete < handle
    % y   : output data
    % u   : input data
    % N   : data length
    % Ts  : sample step
    % N_O : output number
    % N_I : input number
    properties
        y
        u
        N
        Ts
        N_O
        N_I
        option
    end
    
    methods
        function obj = pem_discrete(y,u,Ts)
           obj.y = y;
           obj.u = u;
           obj.N = length(u);
           obj.Ts = Ts;
           obj.N_O = min(size(y));
           obj.N_I = min(size(u));
           % optimiztion option
            obj.option = optimoptions('lsqnonlin');
            obj.option.UseParallel = true;
            obj.option.Display = 'iter-detailed';
            obj.option.FiniteDifferenceType = 'central';
            obj.option.SpecifyObjectiveGradient = true;
            obj.option.StepTolerance = 1e-6;
            obj.option.OptimalityTolerance = 1e-6;
            obj.option.FunctionTolerance = 1e-6;
        end
        
        % params to tf
        function G = get_tf(obj,param)
            G = tf(param,1,obj.Ts,'Variable','q^-1');
        end

        % ARX model (TLS)
        function [a,b] = ARX(obj,dim)
            na = dim(1); nb = dim(2); nk = dim(3);
            Y = zeros(obj.N,na);
            U = zeros(obj.N,nb);
            for i = 1 : na
                Y(:,i) = [zeros(i,1);obj.y(1:end-i)];
            end
            for i = 1 : nb
                U(:,i) = [zeros(i+nk-1,1);obj.u(1:end-nk+1-i)];
            end
            PHI = [-Y(na+1:end,:),U(na+1:end,:)];
            params = PHI\obj.y(na+1:end);
            params = params';
            a = params(1:na);%row
            b = [zeros(1,nk),params(nb:end)];%row
            if nargout < 2
                a = struct(...
                            'A',get_tf(obj,[1,a]),...
                            'B',get_tf(obj,b),...
                            'Ts',obj.Ts...
                            );
            end
        end
        % AR model (TLS)
        function a = AR(obj,na,~)
            condition = (nargin == 2);
            Y = zeros(obj.N,na);
            for i = 1 : na
                Y(:,i) = [zeros(i,1);obj.y(1:end-i)];
            end
            PHI = -Y(na+1:end,:);
            params = PHI\obj.y(na+1:end);
            params = params';
            a = params(1:na);
            if condition
                a = struct(...
                            'A',get_tf(obj,[1,a]),...
                            'Ts',obj.Ts...
                            );
            end
        end
        % ARMAX model (NLS)
        function model = armax(obj,dim,init)
            na = dim(1); nb = dim(2); nc = dim(3);nk = dim(4);
            if nargin == 3
                [b,a] = tfdata(init,'v');
                a = a(2:end);
            else
                [a,b] = obj.ARX([na,nb,nk]);
            end
            params_ini = [a,b,zeros(1,nc)];
            dim = [dim(1:3), 0, 0, dim(end)];
            cost_func = @(theta)obj.pem_discrete_(theta,dim);
            params = lsqnonlin(cost_func,params_ini,[],[],obj.option);
            
            a = params(1:dim(1));
            b = params(dim(1)+1:sum(dim(1:2)));
            c = params(sum(dim(1:2))+1:sum(dim(1:3)));
            model = struct(...
                        'A' , get_tf(obj,[1,a]),...
                        'B' , get_tf(obj,b),...
                        'C' , get_tf(obj,[1,c]),...
                        'Ts', obj.Ts...
                    );
        end
        % OE model (NLS)
        function model = oe(obj,dim,init)
            nb = dim(1); nf = dim(2); nk = dim(3);
            if nargin == 3
                [b,f] = tfdata(init,'v');
                f = f(2:end);
            else
                [f,b] = obj.ARX([nf,nb,nk]);
            end
            params_ini = [b,f];
            dim = [ 0, dim(1), 0, 0, dim(2), dim(3)];
            
            cost_func = @(theta)obj.pem_discrete_(theta,dim);
            params = lsqnonlin(cost_func,params_ini,[],[],obj.option);
            
            b = params(dim(1)+1:sum(dim(1:2)));
            f = params(sum(dim(1:4))+1:sum(dim(1:5)));
            model = struct(...
                        'B' , get_tf(obj,b),...
                        'F' , get_tf(obj,[1,f]),...
                        'Ts', obj.Ts...
                    );
        end
        % PEM model (NLS)
        function model = pem(obj,dim,init)
            na = dim(1); nb = dim(2); nc = dim(3); nd = dim(4); nf = dim(5); nk = dim(6);
            if nargin == 3
                [b,a] = tfdata(init,'v');
                a = a(2:end);
            else
                [a,b] = obj.ARX([na,nb,nk]);
            end
            params_ini = [a,b,zeros(1,nc),zeros(1,nd),zeros(1,nf)];
            
            cost_func = @(theta)obj.pem_discrete_(theta,dim);
            params = lsqnonlin(cost_func,params_ini,[],[],obj.option);
            
            a = params(1:dim(1));
            b = params(dim(1)+1:sum(dim(1:2)));
            c = params(sum(dim(1:2))+1:sum(dim(1:3)));
            d = params(sum(dim(1:3))+1:sum(dim(1:4)));
            f = params(sum(dim(1:4))+1:sum(dim(1:5)));
            model = struct(...
                        'A' , get_tf(obj,[1,a]),...
                        'B' , get_tf(obj,b),...
                        'C' , get_tf(obj,[1,c]),...
                        'D' , get_tf(obj,[1,d]),...
                        'F' , get_tf(obj,[1,f]),...
                        'Ts', obj.Ts...
                    );
        end
        
        % Cost function
        function [f,g] = pem_discrete_(obj,params,dim)
            [yhat,g] = obj.discrete_predictor(params,dim);
            error = obj.y-yhat;
            f = error./obj.N;
        end
        % Predictor
        function [yhat,g] = discrete_predictor(obj,params,dim)
            a = params(1:dim(1));
            b = params(dim(1)+1:sum(dim(1:2)));
            c = params(sum(dim(1:2))+1:sum(dim(1:3)));
            d = params(sum(dim(1:3))+1:sum(dim(1:4)));
            f = params(sum(dim(1:4))+1:sum(dim(1:5)));
            A = obj.get_tf([1,a]);
            B = obj.get_tf(b);
            C = obj.get_tf([1,c]);
            D = obj.get_tf([1,d]);
            F = obj.get_tf([1,f]);
            yhat = lsim( D*B/(C*F), obj.u) + lsim( (1-D*A/C), obj.y);

            g = zeros(obj.N,sum(dim(1:end-1))); % -1 is nk 
            idx = 1;
            % A
            for itr = 1 : dim(1)
                sys_y = -D/C*obj.get_tf([zeros(1,itr),1]);
                g(:,idx) = lsim( sys_y, obj.y);
                idx = idx + 1;
            end
            % B
            for itr = 1 : dim(2)
                sys_u = D/(C*F)*obj.get_tf([zeros(1,itr-1+dim(6)),1]);% B のみ　定数項があるからitr-1
                g(:,idx) = lsim( sys_u, obj.u);
                idx = idx + 1;
            end
            % C
            for itr = 1 : dim(3)
                sys_u = -D*B/(C*C*F)*obj.get_tf([zeros(1,itr+dim(6)),1]);
                sys_y =  D*A/(C*C)*obj.get_tf([zeros(1,itr),1]);
                g(:,idx) = lsim( sys_u, obj.u) +lsim( sys_y, obj.y);
                idx = idx + 1;
            end
            % D
            for itr = 1 : dim(4)
                sys_u =  B/(C*F)*obj.get_tf([zeros(1,itr+dim(6)),1]);
                sys_y = -A/C*obj.get_tf([zeros(1,itr),1]);
                g(:,idx) = lsim( sys_u, obj.u) + lsim( sys_y, obj.y);
                idx = idx + 1;
            end
            % F
            for itr = 1 : dim(5)
                sys_u = -D*B/(C*F*F)*obj.get_tf([zeros(1,itr+dim(6)),1]);
                g(:,idx) = lsim( sys_u, obj.u);
                idx = idx + 1;
            end
        end
        
    end
end