classdef State_Estimator < handle
    % RL state estimator class
    % properties
    % intenal_state : internal state of environmnet state estimator 
    
    properties
        internal_state
    end
    
    methods(Abstract)
        estimate(obj)
        update_internal_state(obj)
    end
end

