classdef Normalized_Radial_Basis_Function < Radial_Basis_Function
    % Normalized Radial Basis Function
    
    properties
        
    end
    
    methods
        function phi = basis_func(obj, x)
            if obj.N > 10000 && obj.N < 2000000% gpuがあると気持ち早くなる．
                phi = exp(-sum((gpuArray(single(x))-gpuArray(single(obj.mu))).^2, 2)./(2*single(obj.sigma).^2));
                phi = phi./sum(phi);
                phi = gather(phi);
            else
                phi = exp(-sum((x-obj.mu).^2, 2))./(2*obj.sigma.^2);
            end
        end
    end
end

