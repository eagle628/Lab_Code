classdef Dynamic_LTI_SS < apx_function
    % Linear time Invariant systm (expressed state space)
    
    properties
        form
        internal_state_c % current step state
        internal_state_p % previous one step state
        grad_internal_state_c
        grad_internal_state_p
    end
    
    methods
        function obj = Dynamic_LTI_SS(gen_ss)
            obj.form = gen_ss;
            obj.theta = obj.form.get_params();
            obj.internal_state_c = zeros(obj.form.n, 1);
            obj.internal_state_p = zeros(obj.form.n, 1);
            obj.grad_internal_state_c = zeros(obj.form.n*(length(obj.theta)), 1);
            obj.grad_internal_state_p = zeros(obj.form.n*(length(obj.theta)), 1);
        end
        
        % override get & set (need for wrapper class)
        function params = get_params(obj)
            params = obj.form.get_params();
%             if ~all(params == obj.theta)
%                 disp('policy dynamic class parameter error');
%             end
        end
        
        function set_params(obj, theta)
            obj.theta = theta;
            obj.form.set_params(theta);
        end
        
        % override superclass initialize
        function initialize(obj)
            obj.internal_state_c = zeros(size(obj.internal_state_c));
            obj.internal_state_p = zeros(size(obj.internal_state_p));
            obj.grad_internal_state_c = zeros(size(obj.grad_internal_state_c));
            obj.grad_internal_state_p = zeros(size(obj.grad_internal_state_p));
        end
        
        function out = predict(obj, state)
            [A, B, C, D] = obj.form.get_ss();
            ABCD = [A,B;C,D];
            xy = ABCD*[obj.internal_state_c;state];
            obj.internal_state_p = obj.internal_state_c;
            obj.internal_state_c  = xy(1 : obj.form.n);
            out = xy(obj.form.n+1 : end);
        end
        
        function grad = grad(obj, state)
            % old version
%             [Abig, Bbig, Cbig, Dbig] = obj.get_sys_big();
%             ABCD = [Abig, Bbig; Cbig, Dbig];
%             xy = ABCD*[obj.internal_state_; obj.grad_internal_state; state'];
%             xx = [obj.internal_state_; obj.grad_internal_state;];
%             x_ = Abig*xx + Bbig*state';
%             y_= Cbig*xx + Dbig*state';
%             xy = [x_; y_];
            [A, ~, C, ~, dA, dB, dC, dD] = obj.form.get_ss();
            dx = reshape(obj.grad_internal_state_c, obj.form.n, []);
            x = obj.internal_state_c;
            Adx = A*dx;
            dAx = vertcat(dA{:})*x;
            Cdx = C*dx;
            dCx = vertcat(dC{:})*x;
            tmpx_ = Adx(:)+dAx+vertcat(dB{:})*state;
            tmpy_ = Cdx(:)+dCx+vertcat(dD{:})*state;
            xy = [tmpx_; tmpy_];
%             tmpx = arrayfun(@(i) double(dA{i})*double(x)+dB{i}*state', 1:numel(dA), 'UniformOutput', false);
%             tmpy = arrayfun(@(i) C*dx(:, i)+dC{i}*x+dD{i}*state', 1:numel(dC), 'UniformOutput', false);
%             xy = vertcat(tmpx{:}, y_);
            obj.grad_internal_state_p = obj.grad_internal_state_c;
            obj.grad_internal_state_c  = xy(1 : length(obj.grad_internal_state_c));
            grad = xy(length(obj.grad_internal_state_c)+1 :  end);
        end
        
        % % Not Used
        function [Abig, Bbig, Cbig, Dbig] = get_sys_big(obj, theta)
            if nargin == 1 || isempty(theta)
                theta = obj.get_params();
            end
            np = length(theta);
            [A, B, C, D, dA, dB, dC, dD] = obj.form.get_ss(theta);
            n = size(A, 1);
            l = size(C, 1);
            % old version
            if false
            Abig = zeros(np*n, (np+1)*n);
            for itr = 1:np
                Abig((itr-1)*n+1:itr*n, itr*n+1:(itr+1)*n) = A;
                if any(any(dA{itr}))
                Abig((itr-1)*n+1:itr*n, 1:n) = dA{itr};
            end
            end
            else
               Abig = []; 
            end
%             Abig = kron(eye(np+1), A);
%             Bbig = kron(ones(np+1, 1), B*0);
%             Bbig(1:n,:) = B;
%             Cbig = kron(eye(np+1), C);
%             Dbig = kron(ones(np+1,1), D*0);
%             Dbig(1:l, :) = D;
%             for itr  = 1:numel(obj.params)
%                 Abig(itr*n+(1:n), 1:n) = dA{itr};
%                 Bbig(itr*n+(1:n), :) = dB{itr};
%                 Cbig(itr*l+(1:l), 1:n) = dC{itr};
%                 Dbig(itr*l+(1:l), :) = dD{itr};
%             end
%             Abig = [cat(1, dA{:}), kron(eye(np), A)];
%             tmp1 = cat(1, dA{:});
%             tmp2 = vertcat(dA{:});
%             tmp3 = kron(eye(np), A);
%             Abig2 = [tmp1, tmp3];
               
%             Abig2 = kron(eye(np), A)
%             for i = 1:numel(dA)
%             
%             end
%             
            Bbig =  cat(1, dB{:});
            Cbig = [cat(1, dC{:}), kron(eye(np), C)];
            Dbig =  cat(1, dD{:});
%             % new version
%             Abig = kron(repmat([1,zeros(1,np)],np+1,1),-A) + kron(eye(np+1), A);
%             Bbig = kron(ones(np+1, 1), B*0);
%             Cbig = kron(repmat([1,zeros(1,np)],np+1,1),-C) + kron(eye(np+1), C);
%             Dbig = kron(ones(np+1,1), D*0);
%             alpha = 1e-2;
%             dA = cellfun(@(x)x*alpha,dA,'UniformOutput',false);
%             dB = cellfun(@(x)x*alpha,dB,'UniformOutput',false);
%             dC = cellfun(@(x)x*alpha,dC,'UniformOutput',false);
%             dD = cellfun(@(x)x*alpha,dD,'UniformOutput',false);
%             Abig = Abig + blkdiag(zeros(size(A)),dA{:});
%             Bbig = Bbig + cat(1, zeros(size(B)),dB{:});
%             Cbig = Cbig + blkdiag(zeros(size(C)),dC{:});
%             Dbig = Dbig + cat(1, zeros(size(D)),dD{:});
%             % grad only
%             Abig = Abig(n+1:end ,:);
%             Bbig = Bbig(n+1:end ,:);
%             Cbig = Cbig(l+1:end ,:);
%             Dbig = Dbig(l+1:end ,:);
        end
    end
end