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
            if obj.N > 10000
                obj.mu = gpuArray(single(mu)); 
            else
                obj.mu = mu;
            end
        end
        
        function phi = basis_func(obj, x)
            if obj.N > 10000 % gpuがあると気持ち早くなる．
                phi = exp(-sum((gpuArray(single(x))-gpuArray(single(obj.mu))).^2, 2)./(2*single(obj.sigma).^2));
                phi = gather(phi);
            else
                phi = exp(-sum((x-obj.mu).^2, 2))./(2*obj.sigma.^2);
            end
        end
        
%         function set_mu(obj, mu)
%             obj.mu = mu;
%         end
%         
%         function set_sigma(obj, sigma)
%             obj.sigma = sigma;
%         end
    end
end

