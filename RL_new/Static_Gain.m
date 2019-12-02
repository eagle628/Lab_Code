classdef Static_Gain < apx_function
		% Static Gain

	properties
	end

	methods
		function obj = Static_Gain(ny, nu, theta)
			if nargin < 3
				theta = zeros(ny, nu);
			end
			obj.theta = theta;
		end

		function out = predict(obj, state)
			out = obj.theta'*state;
		end

		function grad = grad(obj, state)
			grad = state;
		end

		function [a, b, c, d] = get_ss(obj, theta)
			a = 0;
			b = zeros(1, size(obj.theta, 1));
			c = zeros(size(obj.theta, 2), 1);
			d = obj.theta';
		end
	end
end
