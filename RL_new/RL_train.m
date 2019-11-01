classdef RL_train < handle
    % Reinforcement learning train base class
    
    properties(Abstract)
        model
    end
    
    methods(Abstract)
        train(obj)
    end
end

