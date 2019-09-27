classdef TD_lambda < Update_Rule
    % rule of update
    
    properties
        alpha
        lambda
        z
        zeta
    end
    
    methods
        function obj = TD_lambda(alpha, lambda)
            if nargin < 1 || isempty(alpha)
                alpha = 1e-4;
            end
            if nargin < 2 || isempty(lambda)
                lambda = 0.999;
            end
            obj.alpha = alpha;
            obj.lambda = lambda;
            obj.z = 0;
            obj.zeta = 1;
        end
        
        function new_params = updator(obj, params, grad, delta, gamma)
            obj.z = gamma*obj.lambda*obj.z + obj.zeta*grad;
            new_params = params + obj.alpha*delta*obj.z;
            obj.zeta = obj.zeta*gamma;
        end
        
        function initialize(obj)
           obj.zeta = 1;
           obj.z = zeros(size(obj.z));
        end
    end
end

