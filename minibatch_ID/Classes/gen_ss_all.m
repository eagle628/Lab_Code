%% NOTE
% params : ?
% theta  : paramter(double)
% n      : state number
% m      : input number
% l      : output number
%% 
classdef gen_ss_all < gen_ss
    
    properties
        params
        theta
        n
        m
        l
    end
    
    methods
        function obj = gen_ss_all(n, m, l)
            obj.n = n;
            obj.m = m;
            obj.l = l;
            obj.params = cell(3*n-2+n*(m + l) + m*l, 1);
            for itr = 1:n^2
                obj.params{itr} = sprintf('theta_A_%d', itr);
            end
            b = itr;
            for itr = 1:n*m
                obj.params{itr+b} = sprintf('theta_B_%d', itr);
            end
            b = b+itr;
            for itr = 1:l*n
                obj.params{itr+b} = sprintf('theta_C_%d', itr);
            end
            b = b + itr;
            for itr = 1:m*l
                obj.params{itr+b} = sprintf('theta_D_%d', itr);
            end
            obj.theta = zeros(numel(obj.params), 1);
        end
        
        function set_params(obj, theta)
            obj.theta = theta;
        end
        
        function theta = get_params(obj)
            theta = obj.theta;
        end
        
       
        function set_sys(obj, sys)
            Anew = reshape(sys.a, obj.n^2, 1);
            Bnew = reshape(sys.b, obj.n*obj.m, 1);
            Cnew = reshape(sys.c, obj.l*obj.n, 1);
            Dnew = reshape(sys.D, obj.l*lbj.m, 1);
            params_new = [Anew; Bnew; Cnew; Dnew];
            obj.set_params(params_new);
        end
        
        function [A, B, C, D, dA, dB, dC, dD] = get_ss(obj, theta)
            n = obj.n;%#ok
            m = obj.m;%#ok
            l = obj.l;%#ok
            np = numel(obj.params);
            A = reshape(theta(1:n^2), n, n); %#ok
            b = n^2;%#ok
            theta_B = theta(b+(1:n*m));%#ok
            b = b + n*m;%#ok
            theta_C = theta(b+(1:l*n));%#ok
            b = b + n*l;%#ok
            theta_D = theta(b+(1:l*m));%#ok
            B = reshape(theta_B, obj.n, obj.m);
            C = reshape(theta_C, obj.l, obj.n);
            D = reshape(theta_D, obj.l, obj.m);
            if nargout > 4
                dA = cell(np, 1);
                dB = cell(np, 1);
                dC = cell(np, 1);
                dD = cell(np, 1);
                b = 1;
                for itr1 = 1:n%#ok
                    for itr2 = 1:n%#ok
                        M = zeros(n, n);%#ok
                        M(itr1, itr2) = 1;
                        dA{b} = M;
                        dB{b} = B*0;
                        dC{b} = C*0;
                        dD{b} = D*0;
                        b = b+1;
                    end
                end
                b = n^2;%#ok
                for itr = 1:n*m%#ok
                    dA{itr+b} = A*0;
                    M = B*0;
                    M(1+mod(itr-1, n), 1+floor((itr-1)/n)) = 1;%#ok
                    dB{itr+b} = M;
                    dC{itr+b} = C*0;
                    dD{itr+b} = D*0;
                end
                b = b+n*m;%#ok
                for itr = 1:n*l%#ok
                    dA{itr+b} = A*0;
                    dB{itr+b} = B*0;
                    M = zeros(l, n);%#ok
                    M(1+mod(itr-1, l), 1+floor((itr-1)/l)) = 1;%#ok
                    dC{itr+b} = M;
                    dD{itr+b} = D*0;
                end
                b = b+n*l;%#ok
                
                
                for itr = 1:m*l%#ok
                    dA{itr+b} = A*0;
                    dB{itr+b} = B*0;
                    M = zeros(l, m);%#ok
                    M(1+mod(itr-1, l), 1+floor((itr-1)/l)) = 1; %#ok
                    dC{itr+b} = C*0;
                    dD{itr+b} = M;
                end
            end
        end
        
    end
    
end
