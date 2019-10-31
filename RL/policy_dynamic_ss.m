classdef policy_dynamic_ss < policy_class
    % state space ver
    
    properties
        sigma
        internal_state
        internal_state_ % previous one step state
        grad_internal_state
        grad_internal_state_
    end
    
    methods
        function obj = policy_dynamic_ss(gen_ss_canonical, sigma)
            obj.apx_function = gen_ss_canonical;
            obj.sigma = sigma;
            obj.params = obj.apx_function.get_params();
            obj.internal_state  = zeros(obj.apx_function.n, 1);
            obj.internal_state_ = zeros(obj.apx_function.n, 1);
            obj.grad_internal_state  = zeros(obj.apx_function.n*(length(obj.params)), 1);
            obj.grad_internal_state_ = zeros(obj.apx_function.n*(length(obj.params)), 1);
        end
        
        function input = stocastic_policy(obj, state, varargin)
            [a,b,c,d] = obj.apx_function.get_ss(obj.params);
            ABCD = [a,b;c,d];
            xy = ABCD*[obj.internal_state;state'];
            input = xy(size(a, 1)+1 : end);
            obj.internal_state_ = obj.internal_state;
            obj.internal_state  = xy(1 : size(a, 1));
%             tmp = strcmp(varargin, 'Input-Clipping');
%             if sum(tmp)
%                 input(input>=varargin{find(tmp)+1}) = varargin{find(tmp)+1};
%                 input(input<=-varargin{find(tmp)+1}) = -varargin{find(tmp)+1};
%             end
            input = input + obj.sigma^2*randn(size(input));
        end
        
        function input = determistic_policy(obj, state, varargin)
%             [a,b,c,d] = obj.apx_function.get_ss(obj.params);
%             obj.internal_state_ = a*obj.internal_state_ + b*state';
%             input = c*obj.internal_state_ + d*state';
            [a,b,c,d] = obj.apx_function.get_ss(obj.params);
            ABCD = [a,b;c,d];
            xy = ABCD*[obj.internal_state;state'];
            input = xy(size(a, 1)+1 : end);
            obj.internal_state_ = obj.internal_state;
            obj.internal_state  = xy(1 : size(a, 1));
        end
        
        function grad = policy_grad_mu(obj, pre_input, state)
            [~,~,c,d] = obj.apx_function.get_ss(obj.params);
            % old version
%             [Abig, Bbig, Cbig, Dbig] = obj.get_sys_big();
%             ABCD = [Abig, Bbig; Cbig, Dbig];
%             xy = ABCD*[obj.internal_state_; obj.grad_internal_state; state'];
%             xx = [obj.internal_state_; obj.grad_internal_state;];
%             x_ = Abig*xx + Bbig*state';
%             y_= Cbig*xx + Dbig*state';
%             xy = [x_; y_];
            
            [A, B, C, D, dA, dB, dC, dD] = obj.apx_function.get_ss(obj.params);
            dx = reshape(obj.grad_internal_state, size(A, 1), []);
            x = obj.internal_state;
            Adx = A*dx;
            dAx = vertcat(dA{:})*x;
            Cdx = C*dx;
            dCx = vertcat(dC{:})*x;
            tmpx_ = Adx(:)+dAx+vertcat(dB{:})*state';
            tmpy_ = Cdx(:)+dCx+vertcat(dD{:})*state';
            xy = [tmpx_; tmpy_];
%             tmpx = arrayfun(@(i) double(dA{i})*double(x)+dB{i}*state', 1:numel(dA), 'UniformOutput', false);
%             tmpy = arrayfun(@(i) C*dx(:, i)+dC{i}*x+dD{i}*state', 1:numel(dC), 'UniformOutput', false);
%             xy = vertcat(tmpx{:}, y_);
            obj.grad_internal_state_ = obj.grad_internal_state;
            obj.grad_internal_state  = xy(1 : length(obj.grad_internal_state));
            apx_function_grad = xy(length(obj.grad_internal_state)+1 :  end);
            grad = ((pre_input - (c*obj.internal_state_ + d*state'))./(obj.sigma^2))*apx_function_grad;            
        end
        
        function grad = policy_grad_sigma(obj, pre_input, state, theta1, theta2)
            if nargin < 5
                theta2 = obj.get_policy_sigma();
            end
            if nargin < 4 || isempty(theta1)
                theta1 = obj.get_params();
            end
            [~,~,c,d] = obj.apx_function.get_ss(obj.params);
            grad =  ((pre_input - (c*obj.internal_state_ + d*state'))^2 - (obj.sigma).^2)./(obj.sigma) * (1-obj.sigma);
        end
        
        function set_policy_sigma(obj, theta2)
            obj.sigma = theta2;
        end
        
        function sigma = get_policy_sigma(obj)
            sigma = obj.sigma;
        end
        
        function initialize_memory(obj)
            obj.internal_state  = zeros(size(obj.internal_state));
            obj.internal_state_ = zeros(size(obj.internal_state_));
            obj.grad_internal_state  = zeros(size(obj.grad_internal_state));
            obj.grad_internal_state_ = zeros(size(obj.grad_internal_state_));
        end
        
        function [Abig, Bbig, Cbig, Dbig] = get_sys_big(obj, theta)
            if nargin == 1 || isempty(theta)
                theta = obj.get_params();
            end
            np = length(theta);
            [A, B, C, D, dA, dB, dC, dD] = obj.apx_function.get_ss(theta);
            n = size(A, 1);
            l = size(C, 1);
            % old version
            if false
            Abig = zeros(np*n, (np+1)*n);
            for itr = 1:np
                Abig((itr-1)*n+1:itr*n, itr*n+1:(itr+1)*n) = A;
                if any(any(dA{itr}))
                Abig((itr-1)*n+1:itr*n, 1:n) = dA{itr};
            end
            end
            else
               Abig = []; 
            end
%             Abig = kron(eye(np+1), A);
%             Bbig = kron(ones(np+1, 1), B*0);
%             Bbig(1:n,:) = B;
%             Cbig = kron(eye(np+1), C);
%             Dbig = kron(ones(np+1,1), D*0);
%             Dbig(1:l, :) = D;
%             for itr  = 1:numel(obj.params)
%                 Abig(itr*n+(1:n), 1:n) = dA{itr};
%                 Bbig(itr*n+(1:n), :) = dB{itr};
%                 Cbig(itr*l+(1:l), 1:n) = dC{itr};
%                 Dbig(itr*l+(1:l), :) = dD{itr};
%             end
%             Abig = [cat(1, dA{:}), kron(eye(np), A)];
%             tmp1 = cat(1, dA{:});
%             tmp2 = vertcat(dA{:});
%             tmp3 = kron(eye(np), A);
%             Abig2 = [tmp1, tmp3];
               
%             Abig2 = kron(eye(np), A)
%             for i = 1:numel(dA)
%             
%             end
%             
            Bbig =  cat(1, dB{:});
            Cbig = [cat(1, dC{:}), kron(eye(np), C)];
            Dbig =  cat(1, dD{:});
%             % new version
%             Abig = kron(repmat([1,zeros(1,np)],np+1,1),-A) + kron(eye(np+1), A);
%             Bbig = kron(ones(np+1, 1), B*0);
%             Cbig = kron(repmat([1,zeros(1,np)],np+1,1),-C) + kron(eye(np+1), C);
%             Dbig = kron(ones(np+1,1), D*0);
%             alpha = 1e-2;
%             dA = cellfun(@(x)x*alpha,dA,'UniformOutput',false);
%             dB = cellfun(@(x)x*alpha,dB,'UniformOutput',false);
%             dC = cellfun(@(x)x*alpha,dC,'UniformOutput',false);
%             dD = cellfun(@(x)x*alpha,dD,'UniformOutput',false);
%             Abig = Abig + blkdiag(zeros(size(A)),dA{:});
%             Bbig = Bbig + cat(1, zeros(size(B)),dB{:});
%             Cbig = Cbig + blkdiag(zeros(size(C)),dC{:});
%             Dbig = Dbig + cat(1, zeros(size(D)),dD{:});
%             % grad only
%             Abig = Abig(n+1:end ,:);
%             Bbig = Bbig(n+1:end ,:);
%             Cbig = Cbig(l+1:end ,:);
%             Dbig = Dbig(l+1:end ,:);
        end
        
        function update = policy_constraint(obj, new_params, model, belief_N, varargin)
            tmp = strcmp(varargin, 'Invalid-Constraint');
            if sum(tmp)
                update = varargin{find(tmp)+1};
                return
            end
            A = diag(ones(1, belief_N-1), -1);
            A = kron(A, repmat([1,0], model.ny, 1));
            B  = zeros(size(A, 1), model.ny);
            B(1:model.ny,1:model.ny) = eye(model.ny);
            C = A;
            D = B;
            recorder = ss(A,B,C,D,model.Ts);
            target = model.sys_local({'y'},{'u'});
            target = c2d(target, model.Ts);
            [a,b,c,d] = obj.apx_function.get_ss(new_params);
            controller = ss(a,b,c,d,model.Ts);
            true_controller = controller*recorder;
            [Ap,Bp,Cp,~]  = ssdata(target);
            [Ak,Bk,Ck,Dk] = ssdata(true_controller);
            A_all = [Ak,Bk*Cp; Bp*Ck, Ap+Bp*Dk*Cp];
            pole = eig(A_all);
            update = ~(sum(abs(pole)>1));
        end
        
        function update = constraint(obj, new_params, data, varargin)
            model = data.model;
%             AB = data.belief_sys;%varargin{4};
%             A = AB(:, 1:size(AB, 1));
%             B = AB(:, size(AB, 1)+1:end);
%             C = A;
%             D = B;
%             recorder = ss(A,B,C,D,model.Ts);
            recorder = data.recorder;
%             target = model.sys_local({'y','w'},{'u'});
%             target = c2d(target, model.Ts);
            target = model.sys_local_discrete;
            [a,b,c,d] = obj.apx_function.get_ss(new_params);
            [A, B, C, D] = ssdata(recorder);
            Ak = [A, tools.zeros(A, a); b*C, a];
            Bk = [B; b*D];
            Ck = [d*C c];
            Dk = d*D;
%             controller = ss(a,b,c,d,model.Ts);
%             true_controller = controller*recorder;
            [Ap,Bp,Cp,~]  = ssdata(target);
%             [Ak,Bk,Ck,Dk] = ssdata(true_controller);
            A_all = [Ak,Bk*Cp; Bp*Ck, Ap+Bp*Dk*Cp];
            pole = eig(A_all);
            update = abs(max(pole))< 1;
        end
        
        function grad = grad(obj, data)
            grad = obj.policy_grad_mu(data.pre_input, data.state);
        end
    end
end

