classdef gen_ss < handle
    %GEN_SS このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Abstract)
        params
        N
    end
    
    methods(Abstract)
        get_params(obj);
        set_params(obj, theta);
        [A, B, C, D, dA, dB, dC, dD] = get_ss(obj, theta);
    end
    
    methods
        function set_sys(obj, sys)
            
        end
        
        function sys = get_sys(obj, Ts, x)
            if nargin < 3 || isempty(x)
                x = obj.get_params();
            end
           [A, B, C, D] = obj.get_ss(x);
           if nargin < 2 || isempty(Ts)
               sys = ss(A, B, C, D);
           else
               sys = ss(A, B, C, D, Ts);
           end
        end
        function varargout = get_ss_p(obj)
            varargout = cell(nargout, 1);
           [varargout{:}] = obj.get_ss(obj.get_params()); 
        end
        
        function o = order(obj)
            A = obj.get_ss_p();
            o = size(A);
        end
        
        function [c, ceq, dc, dceq] = con_blance(obj, theta)
            [~, B, C, ~, ~, dB, dC] = obj.get_ss(theta);
            c = 0;
            dc = zeros(numel(dB), 1);
            ceq = sum(sum(B.^2)) - sum(sum(C.^2));
            dceq = zeros(numel(dB), 1);
            for itr = 1:numel(dB)
                dceq(itr) = sum(sum(2*dB{itr}.*B))-sum(sum(2*dC{itr}.*C));
            end
        end
        
        function [c, ceq] = con_stable(obj, theta)
            ceq = 0;
            A = obj.get_ss(theta);
            c = max(real(eig(A)));
        end
        
        function [c, ceq, dc, dceq] = constraint(obj, theta)
           c = 0;
           ceq = 0;
           dc = zeros(numel(theta), 1);
           dceq = zeros(numel(theta), 1);
        end
    end
    
end

