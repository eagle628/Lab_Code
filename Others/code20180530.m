% このコードは離れた相手同士の可到達性を確かめるプログラムである．
%% initiallize workspace
clear 
close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
% 20 node 
% node parameter equal
% Cenect ratio 0.1
n = network_swing_simple(30, [1, 2], [2, 10]*1e-2, 1, [1,3], 0.1, 10);
n.Adj_ref = n.Adj_ref*0;
%% local system
%{
com_N = 3;
C = combnk(1:n.N,com_N);
for i = 1 : size(C,1)
    idx = C(i,1:com_N);
    [sys_local, sys_env] = n.get_sys_local(idx);
    A = sys_local.A;
    B = sys_local.B(:,1:numel(idx));
    Co = ctrb(A,B);
    if rank(Co) == numel(idx)*2
        C(i,com_N+1) = 1;
    end
end
%}

idx = [6,15,26,27,28];
[sys_local, sys_env] = n.get_sys_local(idx);
A = sys_local.A;
B = sys_local.B(:,1:numel(idx));
Co = ctrb(A,B);
rank(Co);