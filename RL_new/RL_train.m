classdef RL_train < handle
    % Reinforcement learning train base class
    % model : reinforcement learning environment (inheritanced environment_model class)
    properties
        model
    end
    
    methods(Abstract)
        train(obj)
    end
end

