classdef deep_update_rule < Update_Rule
    %UNTITLED2 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        
    end
    
    methods
        function obj = deep_update_rule()
            
        end
        
        function new_params = updator(obj, params, grad, data, varargin)
            varargin{1}.apx_function.update();
            new_params = [];
        end
    end
end

