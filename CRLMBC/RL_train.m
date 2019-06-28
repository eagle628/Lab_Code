classdef RL_train < handle
    %UNTITLED14 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        model
        policy
        value
    end
    
    methods(Abstract)
        train(obj)
    end
end

