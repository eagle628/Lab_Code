function theta = fit_adam(obj, t, u, y, theta, learning_ratio, rho1, rho2, epsilon, weight)

if nargin < 5 || isempty(theta)
   theta = obj.get_params_fixed();
end

if nargin < 6
    learning_ratio = 1e-3;
end

if nargin < 7
   rho1 = 0.9; 
end

if nargin < 8
    rho2 = 0.999;
end

if nargin < 9
    epsilon = 1e-8;
end

if nargin < 10
    weight = 1;
end

% H = struct();
% H.axes = struct();
% H.f = figure('Name','Callback');
% H.axes.ax1 = subplot(3,1,1);
% H.axes.ax2 = subplot(3,1,2);
% H.axes.ax3 = subplot(3,1,3);

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y(:), @(x) obj.sim_fix(t, u, x));

m = 0;
v = 0;
% n_batch = floor(size(y, 1)/2);
for itr = 1:obj.max_iter
    weight = rand(size(y))>0.8;
%     weight = zeros(size(y));
%     weight(randi(size(y, 1), n_batch), :) = 1;
%     weight(mod(itr, size(y, 1)-1)+1, :) = 1;
    [J, dJ] = obj.eval_func(t, u, y, theta, weight);
    m = rho1*m + (1-rho1)*dJ;
    v = rho2*v + (1-rho2)*(dJ).^2;
    dtheta = -learning_ratio*m/(1-rho1^itr)./sqrt(epsilon+v/(1-rho2^itr));
    theta = theta + dtheta;
    if mod(itr, 10)==0
        J = obj.eval_func(t, u, y, theta);
        func_callback(theta, J, itr);
    end
end

obj.set_params_fixed(theta);

end

