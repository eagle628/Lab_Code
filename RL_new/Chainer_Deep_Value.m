classdef Chainer_Deep_Value < RL_structure
    % apx_function : Chainer optimizer class
    
    properties
        gpu_device
    end
    
    methods
        function obj = Chainer_Deep_Value(apx_function, gpu_device)
            if nargin < 2
                gpu_device = int64(-1);
            end
            obj.apx_function = apx_function;
            obj.gpu_device = gpu_device;
        end
        
        function params = get_params(obj)
          params = [];
        end

        function set_params(obj, params)
          % pass
        end 
        
        function out = predict(obj, state, enable_backprop, varargin)
%             py.chainer.configuration.using_config(pyargs('enable_backprop', enable_backprop))
            py.chainer.using_config.enable_backprop =  enable_backprop;
%             state = py.numpy.array(state,pyargs('dtype','float32'));
%             state = py.numpy.vstack(state);% N*1 2-d numpy array
%             state = state.T;% batch size Single
%             state = py.chainer.dataset.to_device(obj.gpu_device, state);
            state = predict_parser(obj, state);
            out = obj.apx_function.target(state);
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
    end
end
