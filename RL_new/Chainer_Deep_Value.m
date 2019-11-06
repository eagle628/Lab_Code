classdef Chainer_Deep_Value < Value_base
    % apx_function : Chainer optimizer class
    
    properties
        gpu_device
        fixed_apx_function_enable
        fixed_apx_function
    end
    
    methods
        function obj = Chainer_Deep_Value(apx_function, gpu_device, fixed_apx_function_enable)
            if nargin < 3 || isempty(fixed_apx_function_enable)
                fixed_apx_function_enable = false;
            end
            if nargin < 2 || isempty(gpu_device)
                gpu_device = int64(-1);
            end
            obj@Value_base(apx_function)
%             obj.apx_function = apx_function;
            obj.gpu_device = gpu_device;
            obj.fixed_apx_function_enable = fixed_apx_function_enable;
            obj.fixed_apx_function = [];
            if obj.fixed_apx_function_enable
                obj.fixed_apx_function = py.copy.deepcopy(obj.apx_function.target);
            end
        end
        
        function params = get_params(obj)
          params = [];
        end

        function set_params(obj, params)
          % pass
        end 
        
        function out = predict(obj, state, enable_backprop, varargin)
            state = predict_parser(obj, state);
            if obj.fixed_apx_function_enable && ~enable_backprop
                out = obj.fixed_apx_function(state);
                return
            end
            out = obj.apx_function.target.predict(state, enable_backprop);
        end
        
        function grad(obj, data)
            tmp = py.numpy.array(0 ,pyargs('dtype','float32'));
            tmp = tmp.reshape(int64(1),int64(1));
            tmp = py.chainer.dataset.to_device(obj.gpu_device, tmp);
            loss = py.chainer.functions.mean_squared_error(data.delta, tmp);
            obj.apx_function.target.cleargrads();
            loss.backward();
        end
        
        function state = predict_parser(obj, state)
            state = py.numpy.array(state',pyargs('dtype','float32'));
            state = py.numpy.vstack(state);% N*1 2-d numpy array
            state = state.T;% batch size Single
            state = py.chainer.dataset.to_device(obj.gpu_device, state);
        end
        
        function fixed_apx_function_update(obj)
            if obj.fixed_apx_function_enable
                obj.fixed_apx_function = py.copy.deepcopy(obj.apx_function.target);
            end
        end
    end
end
