function [theta ,Jhistory] = fit_adamax(obj, t, u, y, theta, learning_ratio, rho1, rho2, weight, snr)

if nargin < 5 || isempty(theta)
   theta = obj.get_params_fixed();
end

if nargin < 6 || isempty(learning_ratio)
    learning_ratio = 2e-3;
end

if nargin < 7 || isempty(rho1)
   rho1 = 0.9; 
end

if nargin < 8 || isempty(rho2)
    rho2 = 0.999;
end

if nargin < 9 || isempty(weight)
    weight = 1;
end

if nargin < 10
    snr = []; % signal/noise ratio
end

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y(:), @(x) obj.sim_fix(t, u, x));

m = 0;
v = 0;
% n_batch = floor(size(y, 1)/2);
Jhistory = zeros(obj.max_iter, 1);
for itr = 1:obj.max_iter
    weight = rand(size(y))>=(1-weight);
%     weight = zeros(size(y));
%     weight(randi(size(y, 1), n_batch), :) = 1;
%     weight(mod(itr, size(y, 1)-1)+1, :) = 1;
    if ~isempty(snr)
        rng('shuffle')
        d = randn(size(y))*snr;
    else
        d = 0;
    end
    [J, dJ] = obj.eval_func(t, u, y+d, theta, weight);
    Jhistory(itr, 1) = obj.eval_func(t, u, y, theta);
    m = rho1*m + (1-rho1)*dJ;
    v = max(rho2*v, max(abs(dJ)));
    dtheta = -learning_ratio/(1-rho1^itr)*m/v;
    theta = theta + dtheta;
    if mod(itr, 10)==0
        J = obj.eval_func(t, u, y, theta);
        func_callback(theta, J, itr);
    end
%     if itr > 100
%         diff = Jhistory(itr-1) - Jhistory(itr);
%         if (diff > 0) && (diff < 1e-2)
%             break;
%         end
%     end
end

obj.set_params_fixed(theta);

end

