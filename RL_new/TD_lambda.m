classdef TD_lambda < optimizer
    % TD(lambda) update rule
    % lr : learning rate
    % lambda : 
    % z : midterm valirable
    % zeta : dicount factor
    properties
        alpha
        lambda
        z
        zeta
        gamma
    end
    
    methods 
        function obj = TD_lambda(target, alpha, lambda, gamma)
            if nargin < 2 || isempty(alpha)
                alpha = 1e-4;
            end
            if nargin < 3 || isempty(lambda)
                lambda = 0.999;
            end
            if nargin < 4 || isempty(gamma)
                gamma = 1;
            end
            obj@optimizer(target);
            obj.alpha = alpha;
            obj.lambda = lambda;
            obj.z = 0;% 本来は，パラメータの長さ分だけあるが，初回のみブロードキャストされる形式でクリア可能
            obj.zeta = 1;
            obj.gamma = gamma;
        end
        
        function new_params = updator(obj, grad, data)
            if isa(data.delta, 'double')
                delta = data.delta;
            else
                delta = double(py.numpy.squeeze(data.delta.data).tolist);
            end
            obj.z = data.gamma*obj.lambda*obj.z + obj.zeta*grad;
            new_params = obj.target.get_params() + obj.alpha*delta*obj.z;
            obj.zeta = obj.zeta*obj.gamma;
        end
        
        function optimizer_initialize(obj, episode)
           obj.zeta = 1;
           obj.z = zeros(size(obj.z));
           if (episode~=0) && obj.trigger_enable
              if ~mod(episode, obj.trigger_period)
                 obj.alpha = obj.trigger_form(obj.alpha); 
              end
           end
        end
    end
end
