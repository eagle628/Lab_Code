classdef Radial_Basis_Function < apx_function
    % Gausian Radial Basis function Class
    % Arga
    % N : parameter length
    % mu : gausian mean(N*lenght(x))
    % sigma : gausian variance(N*1)
    
    properties
        N
        sigma
        mu
    end
    
    methods
        function obj = Radial_Basis_Function(mu, sigma, theta)
            if ~(length(sigma) == length(mu))
               disp('model parameter length miss ');
            end
            obj.sigma = sigma;
            obj.N = length(sigma);
            if obj.N > 10000 && obj.N < 2000000
                obj.mu = gpuArray(single(mu)); 
            else
                obj.mu = mu;
            end
            if nargin < 3
                theta = zeros(obj.N, 1);
            end
            obj.theta = theta;
        end
        
        function phi = basis_func(obj, x)
            if obj.N > 10000 && obj.N < 2000000% gpuがあると気持ち早くなる．
                phi = exp(-sum((gpuArray(single(x'))-gpuArray(single(obj.mu))).^2, 2)./(2*single(obj.sigma).^2));
                phi = gather(phi);
            else
                phi = exp(-sum((x'-obj.mu).^2, 2))./(2*obj.sigma.^2);
            end
        end
        
        function out = predict(obj, x)
            out = obj.basis_func(x)'*obj.theta;
        end
        
        function grad = grad(obj, state)
            grad = obj.basis_func(state);
        end
        
        function initialize_params(obj)
           obj.theta = zeros(size(obj.theta)); 
        end
    end
end

