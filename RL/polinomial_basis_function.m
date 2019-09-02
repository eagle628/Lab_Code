classdef polinomial_basis_function
    %UNTITLED2 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        N
        mu
        alpha % bias
        beta % dim
    end
    
    methods
        function obj = polinomial_basis_function(N, mu, alpha, beta)
            if nargin < 3 || isempty(alpha)
                alpha = 1;
            end
            if nargin < 4
                beta = 2;
            end
            obj.N = N;
            obj.mu = mu;
            obj.alpha = alpha;
            obj.beta = beta;
        end
        
        function phi = basis_func(obj, x)
            phi = (obj.alpha+sum(x.*obj.mu, 2)).^obj.beta;
        end
    end
end

