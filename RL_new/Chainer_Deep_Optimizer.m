classdef Chainer_Deep_Optimizer < optimizer
    
    properties
    end
    
    methods
        function obj = Chainer_Deep_Optimizer(target)
            obj@optimizer(target);
        end
        
        function opt(obj, data)
            obj.target.grad(data);% deep net backward
            if (~obj.constraint_enable) 
               obj.updator();
               obj.counter = obj.counter + 1;
            end
        end
        
        function updator(obj)
            obj.target.apx_function.update();
        end
        
        % override
        function initialize(obj)
            obj.counter = 0;
            optimizer_initialize(obj);
        end
    end
end