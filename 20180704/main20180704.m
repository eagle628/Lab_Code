%% initiallize workspace
clear 
%close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 4;
seed = 8;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;
n.plot()
%% Set I/O node 
c_n_g = {[3]};
n_n_g = {[2]};
% n_n_g = {[1],[3],[4],[1,4],[1,3],[3,4],[1,3,4]};
%% Signal Power
power = {[1,0]};
%power = {[1,0.1],[1,1],[1,10]};
%% Identification Method
%id_method_g = {'ARMAX','OE'};
id_method_g = {'ARX'};
%% Save directory
location = 'C:\Users\Naoya Inoue\Desktop\Test\2018SummerReport\NEW\noise_less';
%% 
ITR = 10;
%% 
parfor i = 1 : numel(c_n_g)
    c_n = c_n_g{i}
    n_n = n_n_g{1};
    for j = 1 : numel(power)
        id_in_p = power{j}(1);
        noise_p = power{j}(2);
        for k = 1 : numel(id_method_g)
            id_method = id_method_g{k};
            disp('Start______')
            name = func_20180704(n,c_n,n_n,id_in_p,noise_p,id_method,ITR,location);
            disp(strcat('END______',name))
        end
    end
end

beep
