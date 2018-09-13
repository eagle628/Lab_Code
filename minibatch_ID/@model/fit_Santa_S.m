function [theta, Jhistory] = fit_Santa_S(obj, t, u, y, theta, learning_ratio, rho, epsilon, burnin, beta_func, weight_b)
    
    if nargin < 5 || isempty(theta)
       theta = obj.get_params_fixed();
    end

    if nargin < 6 || isempty(learning_ratio)
        learning_ratio = 1e-7;
    end

    if nargin < 7 || isempty(rho)
       rho = 0.999; 
    end

    if nargin < 8 || isempty(epsilon)
        epsilon = 1e-1;
    end

    if nargin < 9 || isempty(burnin)
        burnin = round(obj.max_iter/2);
    end
    
    if nargin < 10 || isempty(beta_func)
        beta_func = @(itr) itr^2;
    end
    
    if nargin < 11 || isempty(weight_b)
        weight_b = 1;
    end
    
    % initialize constant
    C = 100;
    
    % Callback
    func_callback = @(x, v, itr) callback_sgd(x, v, itr, y, @(x) obj.sim_fix(t, u, x));

    % initialize
    rng('shuffle')
    N_params = numel(obj.params);
    mu = sqrt(learning_ratio)*randn(N_params, 1);
    alpha = sqrt(learning_ratio)*C;
    v = zeros(N_params, 1);
    
    Jhistory = zeros(obj.max_iter, 1);
    for itr = 1:obj.max_iter
        beta = beta_func(itr);
        weight = rand(size(y))>=(1-weight_b);% Annealing
        [J, dJ] = obj.eval_func(t, u, y, theta, weight);
        Jhistory(itr, 1) = obj.eval_func(t, u, y, theta);
        % algorithm
        v = rho*v + (1-rho)*(dJ).^2;
        g_c = 1./sqrt(epsilon+sqrt(v));
        theta = theta +g_c.*mu/2;
        if itr == 1
            g_p = g_c;
        end
        if t < burnin
            % exploration
            alpha = alpha + (mu.^2-learning_ratio/beta)./2;
            mu = exp(-alpha/2).*mu ...
                + sqrt(2*g_p*learning_ratio/beta).*randn(N_params,1)...
                + learning_ratio/beta*(1-g_p./g_c)./mu;
            mu = exp(-alpha/2).*mu;
            alpha = alpha + (mu.^22-learning_ratio/beta)./2;
        else
            % refinement
%             alpha = alpha;
            mu = exp(-alpha/2).*mu;
            mu = mu - learning_ratio*(g_c.*dJ);
            mu = exp(-alpha/2).*mu;
        end
        theta = theta + g_c.*mu./2;
        % update
        g_p = g_c;
        if mod(itr, 10)==0
            J = obj.eval_func(t, u, y, theta);
            func_callback(theta, J, itr);
        end
    end
    
obj.set_params_fixed(theta);

end
