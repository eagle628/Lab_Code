%% NOTE
% Only SISO

% params : ?
% theta  : paramter(double)
% n      : state number
% m      : input number
% l      : output number
%% 
classdef gen_ss_canonical < gen_ss
    
    properties
        params
        theta
        n
        m
        l
    end
    
    methods
        function obj = gen_ss_canonical(n, m, l)
            obj.n = n;
            obj.m = m;
            obj.l = l;
            obj.params = cell(n + n + 1, 1);
            for itr = 1:n
                obj.params{itr} = sprintf('theta_A_%d', itr);
            end
            b = itr;
%             for itr = 1:n*(m-1)
%                 obj.params{itr+b} = sprintf('theta_B_%d', itr);
%                 % On canonical form (MIMO), second input B is parameters.
%             end
%             b = b + n*(m-1);
            for itr = 1:l*n
                obj.params{itr+b} = sprintf('theta_C_%d', itr);
                % named by colum direction, Not raw direction. (for "reshape")
            end
            b = b + itr;
            for itr = 1:m*l
                obj.params{itr+b} = sprintf('theta_D_%d', itr);
                % named by colum direction, Not raw direction. (for "reshape")
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
            Anew = sys.a(:,end);
            Bnew = sys.b(:,2:end);
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
            A = [zeros(n-1,1), eye(n-1); reshape(theta(1:n), 1, n)];%#ok
            b = n;%#ok
            if m > 1 %#ok
                theta_B = theta(b+(1:n*(m-1)));%#ok
            else
                theta_B = [];
            end
            b = b + n*(m-1);%#ok
            theta_C = theta(b+(1:l*n));%#ok
            b = b + n*l;%#ok
            theta_D = theta(b+(1:l*m));%#ok
            B = [zeros(n, 1),reshape(theta_B, obj.n, obj.m-1)]; %#ok
            B(end,1) = 1; % Matrix B's [one ,end] is One on  canonnical form.
            C = reshape(theta_C, obj.l, obj.n);
            D = reshape(theta_D, obj.l, obj.m);
            if nargout > 4
                dA = cell(np, 1);
                dB = cell(np, 1);
                dC = cell(np, 1);
                dD = cell(np, 1);
                b = 1;
                for itr = 1:n%#ok
                    M = zeros(n, n);%#ok
                    M(end, itr) = 1;
                    dA{b} = M;
                    dB{b} = B*0;
                    dC{b} = C*0;
                    dD{b} = D*0;
                    b = b+1;
                end
                b = n;%#ok
                for itr = 1:n*(m-1)%#ok
                    dA{itr+b} = A*0;
                    M = B*0;
                    M(1+mod(itr-1, n), 1+1+floor((itr-1)/n)) = 1;%#ok % colum add single not el
                    dB{itr+b} = M;
                    dC{itr+b} = C*0;
                    dD{itr+b} = D*0;
                end
                b = b+n*(m-1);%#ok
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
