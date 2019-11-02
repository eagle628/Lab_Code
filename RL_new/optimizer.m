classdef optimizer < handle
    % optimizer base (optimize for apx_function's parameters)
    % target : RL_structure (Ex. Policy, Value)
    % counter : update number
    
    properties
        target
        counter
    end
    
    methods
        function obj = optimizer(target)
            obj.target = target; 
            obj.counter = 0;
        end
        
        function initialize(obj)
            obj.target.apx_function.initialize();
            obj.counter = 0;
            optimizer_initialize(obj);
        end
        
        function set_up(obj, target)
           obj.target = target; 
        end
        
        function update = constraint(obj, neq_params, data)
            update = true;
            % pass function
        end
        
        function optimizer_initialize(obj)
            % pass function
        end
        
        function opt(obj, data)
            grad = obj.target.grad(data);
            new_params = obj.updator(grad, data);
            if obj.constraint(new_params, data)
               obj.target.set_params(new_params);
               obj.counter = obj.counter + 1;
            end
        end
        
    end
    
    methods(Abstract)
       updator(obj) 
    end
end

