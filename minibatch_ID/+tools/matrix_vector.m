classdef matrix_vector
    properties(Access = private)
        mat
    end
    
    methods
        function obj = matrix_vector(mat)
            obj.mat = mat;
        end
        
        function out = size(obj)
            out = size(obj.mat, 2);
        end
        
        function ind = end(obj, k, n)
            ind = size(obj);
    	end 
        
        function out = subsref(obj, S)
            out = obj.mat(:, S.subs{1});
        end
        
    end
end

