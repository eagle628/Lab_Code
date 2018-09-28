%% NOTE
% params : ?
% theta  : paramter(double)
% n      : state number
% m      : input number
% l      : output number
%% 
classdef gen_ss_tridiag < gen_ss
    
    properties
        params
        theta
        n
        m
        l
    end
    
    methods
        function obj = gen_ss_tridiag(n, m, l)
            obj.n = n;
            obj.m = m;
            obj.l = l;
            obj.params = cell(3*n-2+n*(m + l) + m*l, 1);
            for itr = 1:3*n-2
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
            obj.theta(1:n) = -0.5*(1:n);
            obj.theta(3*n-1:end) = 0.1;
        end
        
        function set_params(obj, theta)
            obj.theta = theta;
        end
        
        function theta = get_params(obj)
            theta = obj.theta;
        end
        
       
        function set_sys(obj, sys)
            A = sys.a;
            Js = get_jordan_blocks(A);
            flag = false;
            for itr = 1:numel(Js)
                J = Js{itr};
                if flag
                    if conj(J(1,1)) ~= J1(1,1)
                       keyboard 
                    end
                    Js{itr} = [];
                    flag = false;
                    continue
                end
                if ~isreal(J(1, 1))
                    d1 = ones(size(J,1)*2-1, 1);
                    d2 = zeros(size(J, 1)*2-1, 1);
                    d1(1:2:end) = imag(J(1,1));
                    d2(1:2:end) = -imag(J(1,1));
                    Js{itr} = diag(ones(size(J, 1)*2, 1)*real(J(1,1))) + diag(d1, 1)+diag(d2, -1);
                    flag = true;
                end
                J1 = J;
            end
            Anew = blkdiag(Js{:});
            t = eye(size(Anew));
            t = t(:);
            options = optimoptions('fsolve', 'MaxFunEvals', 1e4, 'MaxIterations', 2e3,...
                'Display', 'none', 'UseParallel', false, 'Algorithm', 'levenberg-marquardt');
            t0 = fsolve(@(t) [reshape(similar(A, t)-Anew, numel(Anew), 1)], t, options);
            T = reshape(t0, size(Anew));
            Bnew = T\sys.b;
            Cnew = sys.c*T;
%             sys = obj.get_sys;
%             Bnew = sys.b;
%             Cnew = sys.c;
            params_new = [diag(Anew); diag(Anew, -1); diag(Anew, 1); Bnew(:); Cnew(:); sys.d(:)];
            obj.set_params(params_new);
        end
        
        function [A, B, C, D, dA, dB, dC, dD] = get_ss(obj, theta)
            n = obj.n;%#ok
            m = obj.m;%#ok
            l = obj.l;%#ok
            np = numel(obj.params);
            A = diag(theta(1:n))+diag(theta(n+(1:n-1)), -1) + diag(theta(2*n-1+(1:n-1)), 1); %#ok
            b = 3*n-2;%#ok
            theta_B = theta(b+(1:n*m));%#ok
            b = b+ n*m;%#ok
            theta_C = theta(b+(1:l*n));%#ok
            b = b+n*l;%#ok
            theta_D = theta(b+(1:l*m));%#ok
            B = reshape(theta_B, obj.n, obj.m);
            C = reshape(theta_C, obj.l, obj.n);
            D = reshape(theta_D, obj.l, obj.m);
            if nargout > 4
                dA = cell(np, 1);
                dB = cell(np, 1);
                dC = cell(np, 1);
                dD = cell(np, 1);
                for itr = 1:n%#ok
                    M = zeros(n, n);%#ok
                    M(itr, itr) = 1;
                    dA{itr} = M;
                    dB{itr} = B*0;
                    dC{itr} = C*0;
                    dD{itr} = D*0;
                end
                b = n;%#ok
                for itr = 1:n-1%#ok
                    M = zeros(n, n);%#ok
                    M(itr+1, itr) = 1;
                    dA{itr+b} = M;
                    dB{itr+b} = B*0;
                    dC{itr+b} = C*0;
                    dD{itr+b} = D*0;
                end
                b = 2*n-1;%#ok
                for itr = 1:n-1%#ok
                    M = zeros(n, n);%#ok
                    M(itr, itr+1) = 1;
                    dA{itr+b} = M;
                    dB{itr+b} = B*0;
                    dC{itr+b} = C*0;
                    dD{itr+b} = D*0;
                end
                b = 3*n-2;%#ok
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

%% local function
function Js = get_jordan_blocks_old(A)
J = jordan(A);
idx = find(diag(J, 1)~=0);
if isempty(idx)
    Js =  num2cell(diag(J));
elseif numel(idx) == 1
    e = diag(J);
    e(idx:idx+1) = [];
    Js = [{J(idx:idx+1, idx:idx+1)}; num2cell(e)];
else
    e = diag(J);
    d = diff(idx);
    di = find(d~=1);
    di = [0; di; numel(idx)];
    idxes = cell(numel(di)-1, 1);
    for itr = 1:numel(di)-1
        idxes{itr} = [idx(di(itr)+1:di(itr+1)); idx(di(itr+1))+1];
    end
    n = sum(cellfun(@numel, idxes));
    Js = cell(numel(idxes)+size(J, 1)-n, 1);
    for itr = 1:numel(idxes)
        Js{itr} = J(idxes{itr}, idxes{itr});
    end
    e(vertcat(idxes{:})) = [];
    Js(itr+1:end) = num2cell(e);
end
end

function Ahat = similar(A, T)
    T = reshape(T, size(A));
    Ahat = T\(A*T);
end

function Js = get_jordan_blocks(A, epsilon)
if nargin==1
   epsilon = 1e-5; 
end
e = eig(A);
eigs = {};
while ~isempty(e)
    d = abs(e-e(1));
    eigs = [eigs; {e(d<epsilon)}];
    e(d<epsilon) = [];
end
for itr = 1:numel(eigs)
   eigs{itr} = ones(size(eigs{itr}))*mean(eigs{itr}); 
end
Js = cell(numel(eigs), 1);
for itr = 1:numel(eigs)
   Js{itr} = diag(eigs{itr})+diag(ones(numel(eigs{itr})-1, 1), 1); 
end
end
