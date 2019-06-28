classdef value_RBF < approximate_function_class
    %UNTITLED12 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties

    end
    
    methods
        function obj = value_RBF(RBF)
            obj.apx_function = RBF;
            obj.params = zeros(RBF.N, 1);
        end
        
        function value = est_value(obj, state, w)
            obj.set_params(w);
            value = obj.apx_function.basis_func(state)'*w;
        end
        
        function grad = value_grad(obj, state)
           grad =  obj.apx_function.basis_func(state);
        end
    end
end

