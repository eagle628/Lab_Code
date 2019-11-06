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
        function initialize(obj, episode)
            obj.target.apx_function.target.reset_state();% For rnn initialize
            obj.counter = 0;
            optimizer_initialize(obj, episode);
        end
        
        function optimizer_initialize(obj, episode)
           if (episode~=0) && obj.trigger_enable
              if ~mod(episode, obj.trigger_period)
                 obj.target.apx_function.lr = obj.trigger_form(obj.target.apx_function.lr); 
              end
           end
        end
    end
end
