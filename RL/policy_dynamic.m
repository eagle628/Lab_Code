classdef policy_dynamic < policy_class
    %UNTITLED11 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        sigma
        internal_state
        internal_state_ % previous one step state
        grad_internal_state
        grad_internal_state_
    end
    
    methods
        function obj = policy_dynamic(gen_ss_tridiag, sigma)
            obj.apx_function = gen_ss_tridiag;
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
            [~,~,c,d] = obj.apx_function.get_ss(obj.params);
            input = c*obj.internal_state + d*state';
        end
        
        function grad = policy_grad_mu(obj, pre_input, state)
            [~,~,c,d] = obj.apx_function.get_ss(obj.params);
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
            % grad only
            Abig = Abig(n+1:end ,:);
            Bbig = Bbig(n+1:end ,:);
            Cbig = Cbig(l+1:end ,:);
            Dbig = Dbig(l+1:end ,:);
        end
    end
end

