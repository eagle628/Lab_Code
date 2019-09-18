classdef policy_dynamic_tf < policy_class
    % tf ver
    
    properties
        sigma
        internal_state
        internal_state_ % previous one step state
        grad_internal_state
        grad_internal_state_
    end
    
    methods
        function obj = policy_dynamic_tf(tf_class, sigma)
            obj.apx_function = tf_class;
            obj.sigma = sigma;
            obj.params = obj.apx_function.get_params();
            obj.internal_state  = zeros(obj.apx_function.n, 1);
            obj.internal_state_ = zeros(obj.apx_function.n, 1);
            obj.grad_internal_state  = zeros(obj.apx_function.n*(length(obj.params)), 1);
            obj.grad_internal_state_ = zeros(obj.apx_function.n*(length(obj.params)), 1);
        end
        
        function input = stocastic_policy(obj, state, varargin)
            [a,b,c,d] = obj.apx_function.get_sys(obj.params);
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
            [~,~,c,d] = obj.apx_function.get_sys(obj.params);
            input = c*obj.internal_state + d*state';
        end
        
        function grad = policy_grad_mu(obj, pre_input, state)
            [~,~,c,d] = obj.apx_function.get_sys(obj.params);
            % old version
            [Abig, Bbig, Cbig, Dbig] = obj.get_sys_big();
            ABCD = [Abig, Bbig; Cbig, Dbig];
            xy = ABCD*[obj.internal_state_; obj.grad_internal_state; state'];
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
            [~,~,c,d] = obj.apx_function.get_sys(obj.params);
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
%             np = obj.apx_function.N;
            [~, ~, C, D, dA, dB, dC, dD] = obj.apx_function.get_sys(theta);
            n = obj.apx_function.n;
            l = obj.apx_function.l;
            m = obj.apx_function.m;
            Abig21 = cellfun(@(x)x*C, dB(1:n), 'UniformOutput', false);
            Abig = [cat(1, Abig21{:}), blkdiag(dA{1:n})];
            Abig = blkdiag(Abig, blkdiag(dA{n+1:end}));
            Bbig21 = cellfun(@(x)x*D, dB(1:n), 'UniformOutput', false);
            Bbig = cat(1, Bbig21{:},dB{n+1:end});
            Cbig11 = cellfun(@(x)x*C, dD(1:n), 'UniformOutput', false);
            Cbig = [cat(1, Cbig11{:}), blkdiag(dC{1:n})];
            Cbig = blkdiag(Cbig, blkdiag(dC{n+1:end}));
            Dbig11 = cellfun(@(x)x*D, dD(1:n), 'UniformOutput', false);
            Dbig = cat(1, Dbig11{:}, dD{n+1:end});
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
            [a,b,c,d] = obj.apx_function.get_sys(new_params);
            controller = ss(a,b,c,d,model.Ts);
            true_controller = controller*recorder;
            [Ap,Bp,Cp,~]  = ssdata(target);
            [Ak,Bk,Ck,Dk] = ssdata(true_controller);
            A_all = [Ak,Bk*Cp; Bp*Ck, Ap+Bp*Dk*Cp];
            pole = eig(A_all);
            update = ~(sum(abs(pole)>1));
        end
    end
end

