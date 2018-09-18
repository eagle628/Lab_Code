function [theta ,Jhistory] = fit_rmsprop(obj, t, u, y, theta, learning_ratio, rho, epsilon, weight_b, snr)

if nargin < 5 || isempty(theta)
   theta = obj.get_params_fixed();
end

if nargin < 6 || isempty(learning_ratio)
    learning_ratio = 1e-3;
end

if nargin < 7 || isempty(rho)
   rho = 0.99; 
end

if nargin < 8 || isempty(epsilon)
    epsilon = 1e-6;
end

if nargin < 9 || isempty(weight_b)
    weight_b = 1;
end

if nargin < 10 
    snr = []; % signal/noise ratio
end

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y, @(x) obj.sim_fix(t, u, x));

v = 0;
% n_batch = floor(size(y, 1)/2);
Jhistory = zeros(1e4, 1);
rng('shuffle')
for itr = 1:obj.max_iter
    weight = rand(size(y))>=(1-weight_b);
    if ~isempty(snr)
        d = randn(size(y))*snr;
    else
        d = 0;
    end
    [~, dJ] = obj.eval_func(t, u, y+d, theta, weight);
    Jhistory(itr, 1) = obj.eval_func(t, u, y, theta);
    v = rho*v + (1-rho)*(dJ).^2;
    dtheta = -learning_ratio*dJ./sqrt(epsilon+v);
    theta = theta + dtheta;
    if mod(itr, 10)==0
        J = obj.eval_func(t, u, y, theta);
        func_callback(theta, J, itr);
    end
%     if itr > 100
%         diff = Jhistory(itr-1) - Jhistory(itr);
%         if (diff > 0) && (diff < 1)
%             break;
%         end
%     end
%     if mod(itr, 1000) == 0
%         beep
%         prompt = 'Do you want more? Y/N [Y]: ';
%         str = input(prompt,'s');
%         if str == 'n'
%             break;
%         end
%     end

end

obj.set_params_fixed(theta);

end

