classdef approximate_function_class < handle
    %UNTITLED13 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        apx_function
        params
    end
    
    methods
        function set_params(obj, params)
            obj.params = params;
        end
        
        function params = get_params(obj)
            params = obj.params;
        end
        
        function obj = initialize_params(obj)
            obj.params = zeros(size(obj.params));
        end
    end
end

