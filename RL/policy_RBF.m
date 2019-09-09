classdef policy_RBF < policy_class
    %UNTITLED11 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        sigma
    end
    
    methods
        function obj = policy_RBF(RBF, sigma)
            obj.apx_function = RBF;
            obj.sigma = sigma;
            obj.params = zeros(RBF.N, 1);
        end
        
        function input = stocastic_policy(obj, state, theta, varargin)
            if nargin < 3 || isempty(theta)
                theta = obj.get_params();
            else
%                 obj.set_params(theta);
            end
            input = obj.apx_function.basis_func(state)'*theta;
%             input = arrayfun(@obj.apx_function.basis_func, gpuArray(single(state)), gpuArray(single(theta)));
%             input = gather(input);
            tmp = strcmp(varargin, 'Input-Clipping');
            if sum(tmp)
                input(input>=varargin{find(tmp)+1}) = varargin{find(tmp)+1};
                input(input<=-varargin{find(tmp)+1}) = -varargin{find(tmp)+1};
            end
            input = input + obj.sigma^2*randn(size(input));
        end
        
        function input = determistic_policy(obj, state, theta)
            if nargin < 3
                theta = obj.get_params();
            end
            input = obj.apx_function.basis_func(state)'*theta;
        end
        
        function grad = policy_grad_mu(obj, pre_input, state, theta1, theta2)
            if nargin < 5
                theta2 = obj.get_policy_sigma();
            end
            if nargin < 4 || isempty(theta1)
                theta1 = obj.get_params();
            end
            grad = ((pre_input - obj.apx_function.basis_func(state)'*theta1)./(obj.sigma^2))*obj.apx_function.basis_func(state);
%             tmp1 = arrayfun(@obj.determistic_policy, gpuArray(single(state)), gpuArray(single(theta)));
%             tmp2 = arrayfun(@obj.apx_function.basis_func, gpuArray(single(state)), gpuArray(single(theta)));
%             grad = ((pre_input - gather(tmp1))./(obj.sigma)^2)*gather(tmp2);
        end
        
        function grad = policy_grad_sigma(obj, pre_input, state, theta1, theta2)
            if nargin < 5
                theta2 = obj.get_policy_sigma();
            end
            if nargin < 4 || isempty(theta1)
                theta1 = obj.get_params();
            end
            grad =  ((pre_input - obj.apx_function.basis_func(state)'*theta1)^2 - (obj.sigma).^2)./(obj.sigma) * (1-obj.sigma);
        end
        
        function set_policy_sigma(obj, theta2)
%             obj.sigma = 1/(1+exp(-theta2)); % ここをコメントアウトすると確率的政策の分散が初期値に固定される．
        end
        
        function sigma = get_policy_sigma(obj)
            sigma = obj.sigma;
        end
    end
end

