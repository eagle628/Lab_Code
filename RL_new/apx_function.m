classdef apx_function < handle
    % approximate function class
    % theta : paramter
    
    properties
        theta
    end
    
    methods
        function initialize(obj)
           % pass fucntion for override
        end
        
        function params = get_params(obj)
            params = obj.theta;
        end
        
        function set_params(obj, theta)
           obj.theta = theta; 
        end
    end
    
    methods(Abstract)
       predict(obj)
       grad(obj)
    end
end

