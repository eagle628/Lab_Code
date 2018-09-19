function [theta ,Jhistory] = fit_adadelta(obj, t, u, y, theta, rho, epsilon, weight_b)

if nargin < 5 || isempty(theta)
   theta = obj.get_params_fixed();
end

if nargin < 6 || isempty(rho)
   rho = 0.95; 
end

if nargin < 7 || isempty(epsilon)
    epsilon = 1e-8;
end

if nargin < 8 || isempty(weight_b)
    weight_b = 1;
end

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y, @(x) obj.sim_fix(t, u, x));

mu = 0;
v = 0;
% n_batch = floor(size(y, 1)/2);
Jhistory = zeros(1e4, 1);
rng('shuffle')
for itr = 1:obj.max_iter
    weight = rand(size(y))>=(1-weight_b);
    [~, dJ] = obj.eval_func(t, u, y, theta, weight);
    Jhistory(itr, 1) = obj.eval_func(t, u, y, theta);
    v = rho*v + (1-rho)*dJ.^2;
    dtheta = -sqrt(mu+epsilon).*dJ./sqrt(v+epsilon);
    mu = rho*mu + (1-rho)*dtheta.^2;
    theta = theta + dtheta;
    if mod(itr, 10)==0
        J = obj.eval_func(t, u, y, theta);
        func_callback(theta, J, itr);
    end
end

obj.set_params_fixed(theta);

end

