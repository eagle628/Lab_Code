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
        
        function input = stocastic_policy(obj, state, theta)
            obj.set_params(theta);
            input = obj.apx_function.basis_func(state)'*theta;
            input = input + obj.sigma*randn(size(input));
        end
        
        function input = determistic_policy(obj, state, theta)
            input = obj.apx_function.basis_func(state)'*theta;
        end
        
        function grad = policy_grad(obj, pre_input, state, theta)
             grad = ((pre_input - obj.apx_function.basis_func(state)'*theta)./(obj.sigma)^2)*obj.apx_function.basis_func(state);
        end
        
    end
end

