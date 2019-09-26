clear
close all

n = 2;
l = 1;
m = 1;

sys = gen_ss_canonical(n,m,l);
init_sys = drss(n,l,m);
sys.set_sys(init_sys);

N = 1000;
Ts = 0.01;
t = (0:N-1)'*Ts;
ddd = randn(N, m);

[a,b,c,d,da,db,dc,dd] = sys.get_ss(sys.theta);

yyy = lsim(init_sys, ddd);
idx = 1;
diff_sys1 = cell(sys.N, 1);
diff_yyy1 = cell(sys.N, 1);
for k = 1 : sys.n
    shape_c = zeros(size(c));
    shape_c(k) = -1;
    ABCD = [a,b;c,zeros(size(d))] * [a,b;shape_c,zeros(size(d))];
    diff_sys1{idx} = ss(ABCD(1:sys.n,1:sys.n),ABCD(1:sys.n, sys.n+1:end),...
                ABCD(sys.n+1:end,1:sys.n),ABCD(sys.n+1:end, sys.n+1:end), Ts);
    diff_yyy1{idx} = lsim(diff_sys1{idx}, ddd);
    idx = idx + 1;
end

for k = 1 : sys.l+1
    shape_c = zeros(size(c));
    shape_c(k) = 1;
    ABCD = [a,b;shape_c,zeros(size(d))];
    diff_sys1{idx} = ss(ABCD(1:sys.n,1:sys.n),ABCD(1:sys.n, sys.n+1:end),...
                ABCD(sys.n+1:end,1:sys.n),ABCD(sys.n+1:end, sys.n+1:end), Ts);
    diff_yyy1{idx} = lsim(diff_sys1{idx}, ddd);
    idx = idx + 1;
end

for k = 1 : sys.m
    ABCD = [zeros(size(a)),zeros(size(b));zeros(size(c)),1];
    diff_sys1{idx} = ss(ABCD(1:sys.n,1:sys.n),ABCD(1:sys.n, sys.n+1:end),...
                ABCD(sys.n+1:end,1:sys.n),ABCD(sys.n+1:end, sys.n+1:end), Ts);
    diff_yyy1{idx} = lsim(diff_sys1{idx}, ddd);
    idx = idx + 1;
end

figure
hold on
cellfun(@plot, repmat({t},sys.N,1), diff_yyy1)

[Abig, Bbig, Cbig, Dbig] = get_sys_big(sys);
diff_sys2 = ss(Abig,Bbig,Cbig,Dbig,Ts);
diff_yyy2 = lsim(diff_sys2, ddd);
diff_yyy2 = mat2cell(diff_yyy2', ones(1,sys.N+1), N);

figure
hold on
cellfun(@plot, repmat({t},sys.N,1), diff_yyy2(2:end))
%% local

function [Abig, Bbig, Cbig, Dbig] = get_sys_big(gen_ss, theta)
            if nargin == 1 || isempty(theta)
                theta = gen_ss.get_params();
            end
            np = length(theta);
            [A, B, C, D, dA, dB, dC, dD] = gen_ss.get_ss(theta);
            n = size(A, 1);
            l = size(C, 1);
%             % old version
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
            % new version
            Abig = kron(repmat([1,zeros(1,np)],np+1,1),-A) + kron(eye(np+1), A);
            Abig(1:n, 1:n) = A;
            Bbig = kron(ones(np+1, 1), B*0);
            Bbig(1:n, :) = B;
            Cbig = kron(repmat([1,zeros(1,np)],np+1,1),-C) + kron(eye(np+1), C);
            Cbig(1:l, 1:n) = C;
            Dbig = kron(ones(np+1,1), D*0);
            Dbig(1:l, :) = D;
            alpha = 1e-3;
            dA = cellfun(@(x)x*alpha,dA,'UniformOutput',false);
            dB = cellfun(@(x)x*alpha,dB,'UniformOutput',false);
            dC = cellfun(@(x)x*alpha,dC,'UniformOutput',false);
            dD = cellfun(@(x)x*alpha,dD,'UniformOutput',false);
            Abig = Abig + blkdiag(zeros(size(A)),dA{:});
            Bbig = Bbig + cat(1, zeros(size(B)),dB{:});
            Cbig = Cbig + blkdiag(zeros(size(C)),dC{:});
            Dbig = Dbig + cat(1, zeros(size(D)),dD{:});
%             % grad only
%             Abig = Abig(n+1:end ,:);
%             Bbig = Bbig(n+1:end ,:);
%             Cbig = Cbig(l+1:end ,:);
%             Dbig = Dbig(l+1:end ,:);
        end
