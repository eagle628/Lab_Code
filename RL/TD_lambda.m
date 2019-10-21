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
            obj.z = 0;% 本来は，パラメータの長さ分だけあるが，初回のみブロードキャストされる形式でクリア可能
            obj.zeta = 1;
        end
        
        function new_params = updator(obj, params, grad, data, varargin)
            obj.z = data.gamma*obj.lambda*obj.z + obj.zeta*grad;
            if strcmp(class(data.delta), 'couble')
                delta = data.delta;
            else
                delta = single(py.numpy.squeeze(data.delta).tolist);
            end
            new_params = params + obj.alpha*delta*obj.z;
            obj.zeta = obj.zeta*data.gamma;
        end
        
        function initialize(obj)
           obj.zeta = 1;
           obj.z = zeros(size(obj.z));
        end
    end
end

