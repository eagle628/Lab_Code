function [theta ,Jhistory] = fit_adasecant(obj, t, u, y, theta, learning_ratio)

if nargin < 5 || isempty(theta)
   theta = obj.get_params_fixed();
end

if nargin < 6 || isempty(learning_ratio)
    learning_ratio = 1e-3;
end


end

