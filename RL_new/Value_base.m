classdef Value_base < RL_structure
   % Stocastic policy class
   % Args
   % apx_funciton : apx_funciton class
   
   properties
      
   end
   
   methods
       function obj = Value_base(apx_function)
          obj.apx_function = apx_function;
       end
       
       function out = predict(obj, state, varargin)
           out = obj.apx_function.predict(state);
       end
       
       function grad = grad(obj, data)
           grad = obj.apx_function.grad(data.state);
       end
   end
end
