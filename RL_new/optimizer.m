classdef optimizer < handle
    % optimizer base (optimize for apx_function's parameters)
    % target : RL_structure (Ex. Policy, Value)
    % counter : update number
    
    properties
        target
        constraint_enable
        counter
        trigger_enable
        trigger_period
        trigger_form
    end
    
    methods
        function obj = optimizer(target, constraint_enable)
            if nargin < 2
                constraint_enable = false;
            end
            obj.target = target; 
            obj.counter = 0;
            obj.constraint_enable = constraint_enable;
            obj.trigger_enable = false;
            obj.trigger_period = 100;
            obj.trigger_form = @(x)0.1*x;
        end
        
        function initialize(obj, episode)
            obj.target.apx_function.initialize();
            obj.counter = 0;
            optimizer_initialize(obj, episode);
        end
        
        function set_up(obj, target)
           obj.target = target; 
        end
        
        function update = constraint(obj, new_params, data)
            update = data.model.constraint(obj.target, new_params, data);
        end
        
        function optimizer_initialize(obj, episode)
            % pass function
        end
        
        function grad = opt(obj, data)
            grad = obj.target.grad(data);
            new_params = obj.updator(grad, data);
            if obj.constraint_enable
                if obj.constraint(new_params, data)
                    grad = double(py.numpy.squeeze(data.delta.data).tolist)*grad;
                    obj.target.set_params(new_params);
                    obj.counter = obj.counter + 1;
                else
                    grad = grad*0;
                end
            else
                grad = data.delta*grad;
                obj.target.set_params(new_params);
                obj.counter = obj.counter + 1;
            end
        end
        
    end
    
    methods(Abstract)
       updator(obj) 
    end
end

