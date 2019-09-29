classdef value_RBF < value_class
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
            if nargin < 3 || isempty(w)
                w = obj.get_params();
            else
%                 obj.set_params(w);
            end
            value = obj.apx_function.basis_func(state)'*w;
%             value = arrayfun(@obj.apx_function.basis_func, gpuArray(single(state)), gpuArray(single(w)));
%             value = gather(value);
        end
        
        function grad = value_grad(obj, state, varargin)
            grad =  obj.apx_function.basis_func(state);
            tmp = strcmp(varargin, 'Grad-Clipping');
            if sum(tmp)
                grad(grad>=varargin{find(tmp)+1}) = varargin{find(tmp)+1};
                grad(grad<=-varargin{find(tmp)+1}) = -varargin{find(tmp)+1};
            end
        end
        
        function grad = grad(obj, state, varargin)
           grad = obj.value_grad(state, varargin{:}); 
        end
        
        function value = predict(obj, state, w)
            if nargin < 3
                w = obj.get_params();
            end
            value = obj.est_value(state, w); 
        end
    end
end

