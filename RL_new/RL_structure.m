classdef RL_structure < handle
    % Base class for value or policy
    % Arags
    % optimizer : approximate function
    
    properties
        apx_function
    end
    
    methods(Abstract)
       grad(obj)
       predict(obj)
    end
    
    methods
       function params = get_params(obj)
            params = obj.apx_function.get_params();
        end
        
        function set_params(obj, params)
            obj.apx_function.set_params(params);
        end
    end
end

