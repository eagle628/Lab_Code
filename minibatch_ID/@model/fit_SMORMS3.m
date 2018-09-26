function [theta ,Jhistory] = fit_SMORMS3(obj, t, u, y, theta, learning_ratio, epsilon, weight_b)

if nargin < 5 || isempty(theta)
   theta = obj.get_params_fixed();
end

if nargin < 6 || isempty(learning_ratio)
    learning_ratio = 1e-3;
end

if nargin < 7 || isempty(epsilon)
    epsilon = 1e-6;
end

if nargin < 8 || isempty(weight_b)
    weight_b = 1;
end

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y, @(x) obj.sim_fix(t, u, x));

% n_batch = floor(size(y, 1)/2);
Jhistory = zeros(1e4, 1);
rng('shuffle')
s = 0;
m = 0;
v = 0;
zeta = 0;
learning_ratio = ones(size(theta))*learning_ratio;
for itr = 1:obj.max_iter
    weight = rand(size(y))>=(1-weight_b);
    [~, dJ] = obj.eval_func(t, u, y, theta, weight);
    Jhistory(itr, 1) = obj.eval_func(t, u, y, theta);
    s = 1 + (1-zeta.*s);
    rho = 1./(s+1);
    m = rho.*m + (1-rho).*dJ;
    v = rho.*v + (1-rho).*(dJ).^2;
    zeta = m.^2./(v+epsilon);
    dtheta = - min([learning_ratio,zeta],[],2)./(v+epsilon).*dJ;
    theta = theta + dtheta;
    if mod(itr, 10)==0
        J = obj.eval_func(t, u, y, theta);
        func_callback(theta, J, itr);
    end
end

obj.set_params_fixed(theta);

end

