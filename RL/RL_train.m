classdef RL_train < handle
    %UNTITLED14 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        model
    end
    
    methods(Abstract)
        train(obj)
        reward(obj)
    end
end

