classdef optimizer < handle
    %UNTITLED3 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        update_rule
        approximate_function_class
        counter
    end
    
    methods
        function obj = optimizer(update_rule, approximate_function_class)
            obj.update_rule = update_rule;
            obj.approximate_function_class = approximate_function_class;
            obj.counter = 0;
        end
        
        function opt(obj, data)
            grad = obj.approximate_function_class.grad(data);
            params = obj.approximate_function_class.get_params();
            new_params = obj.update_rule.updator(params, grad, data, obj.approximate_function_class);
            if obj.approximate_function_class.constraint(new_params, data)
                obj.approximate_function_class.set_params(new_params);
                obj.counter = obj.counter + 1;
            end
        end
        
        function initialize(obj)
            obj.approximate_function_class.initialize_memory();
            obj.update_rule.initialize();
            obj.counter = 0;
        end
    end
end

