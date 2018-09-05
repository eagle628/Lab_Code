function [outputArg1,outputArg2] = fit_momentum(obj, t, u, y, theta, learning_ratio, mu, weight)

if nargin < 5
   theta = obj.get_params_fixed();
end

if nargin < 6
    learning_ratio = 1e-4;
end

if nargin < 7
   mu = 0.95; 
end

if nargin < 8
    weight = 1;
end

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y(:), @(x) obj.sim_fix(t, u, x));

dtheta= 0;
for itr = 1:obj.max_iter
%     weight = rand(size(y));
    [v, dv] = obj.eval_func(t, u, y, theta, weight);
    dtheta = mu*dtheta - (1-mu)*learning_ratio*dv;
    theta = theta + dtheta;
%     if mod(itr, 100)==1
        func_callback(theta, v, itr);
%     end
end

obj.set_params_fixed(theta);

end

