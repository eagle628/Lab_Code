classdef policy_RBF < approximate_function_class
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
            obj.set_params(theta);
            input = obj.apx_function.basis_func(state)'*theta;
%             input = arrayfun(@obj.apx_function.basis_func, gpuArray(single(state)), gpuArray(single(theta)));
%             input = gather(input);
            tmp = strcmp(varargin, 'clipping');
            if sum(tmp)
                if  strcmp(varargin{find(tmp)+1},'on')
                    input(input>=10) = 10;
                    input(input<=-10) = -10;
                end
            end
            input = input + obj.sigma^2*randn(size(input));
        end
        
        function input = determistic_policy(obj, state, theta)
            input = obj.apx_function.basis_func(state)'*theta;
        end
        
        function grad = policy_grad_mu(obj, pre_input, state, theta1, theta2)
            narginchk(4, inf)
            grad = ((pre_input - obj.apx_function.basis_func(state)'*theta1)./(obj.sigma^2))*obj.apx_function.basis_func(state);
%             tmp1 = arrayfun(@obj.determistic_policy, gpuArray(single(state)), gpuArray(single(theta)));
%             tmp2 = arrayfun(@obj.apx_function.basis_func, gpuArray(single(state)), gpuArray(single(theta)));
%             grad = ((pre_input - gather(tmp1))./(obj.sigma)^2)*gather(tmp2);
        end
        
        function grad = policy_grad_sigma(obj, pre_input, state, theta1, theta2)
            narginchk(4, inf)
            grad =  ((pre_input - obj.apx_function.basis_func(state)'*theta1)^2 - (obj.sigma).^2)./(obj.sigma) * (1-obj.sigma);
        end
        
        function set_policy_sigma(obj, theta2)
            obj.sigma = 1/(1+exp(-theta2)); % ここをコメントアウトすると確率的政策の分散が初期値に固定される．
        end
        
        function sigma = get_policy_sigma(obj)
            sigma = obj.sigma;
        end
    end
end

