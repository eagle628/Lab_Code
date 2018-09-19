function [theta ,Jhistory] = fit_adam(obj, t, u, y, theta, learning_ratio, rho1, rho2, epsilon, weight_b)

if nargin < 5 || isempty(theta)
   theta = obj.get_params_fixed();
end

if nargin < 6 || isempty(learning_ratio)
    learning_ratio = 1e-3;
end

if nargin < 7 || isempty(rho1)
   rho1 = 0.9; 
end

if nargin < 8 || isempty(rho2)
    rho2 = 0.999;
end

if nargin < 9 || isempty(epsilon)
    epsilon = 1e-8;
end

if nargin < 10 || isempty(weight_b)
    weight_b = 1;
end

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y, @(x) obj.sim_fix(t, u, x));

m = 0;
v = 0;
% n_batch = floor(size(y, 1)/2);
Jhistory = zeros(1e4, 1);
rng('shuffle')
for itr = 1:obj.max_iter
    weight = rand(size(y))>=(1-weight_b);
%     weight = zeros(size(y));
%     weight(randi(size(y, 1), n_batch), :) = 1;
%     weight(mod(itr, size(y, 1)-1)+1, :) = 1;
    [~, dJ] = obj.eval_func(t, u, y, theta, weight);
    Jhistory(itr, 1) = obj.eval_func(t, u, y, theta);
    m = rho1*m + (1-rho1)*dJ;
    v = rho2*v + (1-rho2)*(dJ).^2;
    dtheta = -learning_ratio*m/(1-rho1^itr)./sqrt(epsilon+v/(1-rho2^itr));
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

