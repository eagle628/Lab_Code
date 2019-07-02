classdef Radial_Basis_Function
    %UNTITLED10 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        N
        sigma
        mu
    end
    
    methods
        function obj = Radial_Basis_Function(N, mu, sigma)
            obj.sigma = sigma;
            obj.N = N;
            obj.mu = mu; 
        end
        
        function phi = basis_func(obj, x)
           phi = exp(-sum((gpuArray(single(x))-gpuArray(single(obj.mu))).^2, 2)./(2*single(obj.sigma).^2));
           phi = gather(phi);
        end
    end
end

