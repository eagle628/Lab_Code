classdef value_chainer_deep_net < value_class
    %UNTITLED このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        
    end
    
    methods
        function obj = value_chainer_deep_net(deep_net_opt)
            obj.apx_function = deep_net_opt;%  chainer optimizer class
            obj.params = cell(1);
        end
        
        function grad = value_grad(obj, data, varargin)
            tmp = py.numpy.array(0 ,pyargs('dtype','float32'));
%             tmp = py.numpy.vstack(tmp);
%             tmp = tmp.T;
            tmp = tmp.reshape(int64(1),int64(1));
            tmp = py.chainer.dataset.to_device(int64(0), tmp);
            loss = py.chainer.functions.mean_squared_error(data.delta, tmp);
            obj.apx_function.target.cleargrads();
            loss.backward();
            grad = [];
        end
        
        function value = est_value(obj, state, w)
            if nargin < 3 ||isempty(w)
                w = obj.get_params();
            end
%             value = obj.apx_function.target(py.cupy.array(state,pyargs('dtype','float32')));
            tmp = py.numpy.array(state,pyargs('dtype','float32'));
            tmp = py.numpy.vstack(tmp);
            tmp = tmp.T;
%             tmp.to_device(int64(0))
            tmp = py.chainer.dataset.to_device(int64(0), tmp);
            value = obj.apx_function.target(tmp);
            value = value.data;
        end
        
        function grad = grad(obj, data, varargin)
            grad = obj.value_grad(data, varargin);
        end
    end
end

