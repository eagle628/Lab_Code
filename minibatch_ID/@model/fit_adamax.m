function [theta ,Jhistory] = fit_adamax(obj, t, u, y, theta, learning_ratio, rho1, rho2, weight)

if nargin < 5 || isempty(theta)
   theta = obj.get_params_fixed();
end

if nargin < 6
    learning_ratio = 2e-3;
end

if nargin < 7
   rho1 = 0.9; 
end

if nargin < 8
    rho2 = 0.999;
end

if nargin < 9
    weight = 1;
end

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y(:), @(x) obj.sim_fix(t, u, x));

m = 0;
v = 0;
% n_batch = floor(size(y, 1)/2);
Jhistory = zeros(obj.max_iter, 1);
for itr = 1:obj.max_iter
    weight = rand(size(y))>0.8;
%     weight = zeros(size(y));
%     weight(randi(size(y, 1), n_batch), :) = 1;
%     weight(mod(itr, size(y, 1)-1)+1, :) = 1;
    [J, dJ] = obj.eval_func(t, u, y, theta, weight);
    Jhistory(itr, 1) = J;
    m = rho1*m + (1-rho1)*dJ;
    v = max(rho2*v, abs(dJ));
    dtheta = -learning_ratio/(1-rho1^itr)*m/v;
    theta = theta + dtheta;
    if mod(itr, 10)==0
        J = obj.eval_func(t, u, y, theta);
        func_callback(theta, J, itr);
    end
    if itr > 100
        diff = Jhistory(itr-1) - Jhistory(itr);
        if (diff > 0) && (diff < 1e-4)
            break;
        end
    end
end

obj.set_params_fixed(theta);

end

