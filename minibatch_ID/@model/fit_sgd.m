function [outputArg1,outputArg2] = fit_sgd(obj, t, u, y, theta, learning_ratio, weight)

if nargin < 5
   theta = obj.get_params_fixed();
end

if nargin < 6
    learning_ratio = 1e-4;
end

if nargin < 7
    weight = 1;
end

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y(:), @(x) obj.sim_fix(t, u, x));


for itr = 1:obj.max_iter
%     weight = rand(size(y));
    [v, dv] = obj.eval_func(t, u, y, theta, weight);
    theta = theta - learning_ratio * dv;
%     if mod(itr, 100)==1
        func_callback(theta, v, itr);
%     end
end

obj.set_params_fixed(theta);

end

