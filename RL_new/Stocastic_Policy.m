classdef Stocastic_Policy < RL_structure
   % Stocastic policy class
   % Args
   % apx_funciton : apx_funciton class
   % pi_sigma : stcastic policy variance
   
   properties
      pi_sigma 
      pi_grad_enable
   end
   
   methods
        function obj = Stocastic_Policy(apx_function, pi_sigma)
            obj.apx_function = apx_function;
            obj.pi_sigma = 1/(1+exp(-pi_sigma));
            obj.pi_grad_enable = false;
        end
       
        function [input_real, input] = predict(obj, state, random)
            if nargin < 3
                random = true;
            end
            input = obj.apx_function.predict(state);
            if random
                input_real = input + obj.pi_sigma^2*randn(size(input));
            else
                input_real = input;
            end
        end
       
        function grad = grad(obj, data)
            grad = ((data.pre_input - data.pre_input_mu)./(obj.pi_sigma^2))*obj.apx_function.grad(data.state);
            if obj.pi_grad_enable
                grad2 = obj.grad2(data);
            else
                grad2 = 0;
            end
            grad = [grad; grad2];
        end
        
        function grad = grad2(obj, data)
            grad =  ((data.pre_input - data.pre_input_mu)^2 - (obj.pi_sigma).^2)./(obj.pi_sigma) * (1-obj.pi_sigma);
        end
       
        function set_policy_sigma(obj, theta2)
            obj.pi_sigma = 1/(1+exp(-theta2)); 
        end
       
        function sigma = get_policy_sigma(obj)
            sigma = obj.pi_sigma;
        end
        
        function params = get_params(obj)
           params = [obj.apx_function.get_params(); obj.get_policy_sigma];
        end
        
        function set_params(obj, params)
           obj.apx_function.set_params(params(1:end-1)) ;
           obj.set_policy_sigma(params(end));
        end
    end
end
