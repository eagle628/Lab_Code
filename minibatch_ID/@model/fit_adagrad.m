function [theta ,Jhistory] = fit_adagrad(obj, t, u, y, theta, learning_ratio, weight_b)

if nargin < 5 || isempty(theta)
   theta = obj.get_params_fixed();
end

if nargin < 6 || isempty(learning_ratio)
    learning_ratio = 1e-3;
end

if nargin < 7 || isempty(weight_b)
    weight_b = 1;
end

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y, @(x) obj.sim_fix(t, u, x));

s = 0; % sum of gradient vector
Jhistory = zeros(1e4, 1);
rng('shuffle')
for itr = 1:obj.max_iter
    weight = rand(size(y))>=(1-weight_b);
    [~, dJ] = obj.eval_func(t, u, y, theta, weight);
    Jhistory(itr, 1) = obj.eval_func(t, u, y, theta);
    s = dJ.^2;
    dtheta = -learning_ratio*dJ./sqrt(s);
    theta = theta + dtheta;
    if mod(itr, 10)==0
        J = obj.eval_func(t, u, y, theta);
        func_callback(theta, J, itr);
    end
end

obj.set_params_fixed(theta);

end

