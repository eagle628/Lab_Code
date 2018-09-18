%% Note
% Levenberg–Marquardt algorithm
%% Main
function [theta, Jhistory] = fit_LMA(obj, t, u, y, theta, lambda, weight)
if nargin < 5 || isempty(theta)
   theta = obj.get_params_fixed();
end
if nargin < 6 || isempty(lambda)
   lambda = 1;
end
if nargin < 7 || isempty(weight)
   weight = ones(size(y));
end

func_callback = @(x, v, itr) callback_sgd(x, v, itr, y, @(x) obj.sim_fix(t, u, x));

Jhistory = zeros(obj.max_iter, 1);
for itr = 1:obj.max_iter
    [J, dJ] = obj.eval_func(t, u, y, theta, weight);
%     Jhistory(itr, 1) = obj.eval_func(t, u, y, theta);
    ddJ = diag(dJ)*kron(ones(size(dJ)),dJ');
    while true
        alpha = ddJ + diag(lambda*dJ.^2);
        beta = J.*dJ;
        dtheta = alpha\beta;
        theta_p = theta + dtheta;
        % update lambda
        J_p = obj.eval_func(t, u, y, theta_p);
        if J_p < J
            lambda = lambda*0.5;
            theta = theta_p;
            break;
        end
        lambda = lambda*2;
    end
    % callback
    J = obj.eval_func(t, u, y, theta);
    func_callback(theta, J, itr);
end

obj.set_params_fixed(theta);


end
