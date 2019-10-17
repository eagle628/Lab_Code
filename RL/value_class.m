classdef value_class < approximate_function_class
    %UNTITLED7 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
    end
    
    methods(Abstract)
        est_value(obj)
        grad(obj)
    end
end

