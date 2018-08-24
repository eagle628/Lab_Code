close all
clear
    
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)


%% Set I/O node 
c_n = 1;
n_n = [2,3,4];
%% Signal Power
id_in_p = 1;
noise_p = 1;
%%
% seed_g = [2,4,8,10];
seed_g = 4;
%% save directory
location ='C:\Users\Naoya Inoue\Desktop\Test\2018SummerReport\SPEM';
%% Iteration
ITR = 1;
sim_seed = randi(10000000, ITR, 1);


for i = 1
    if i == 1 
        id_method = 'ARX';
    elseif i == 2
        id_method = 'ARMAX';
    elseif i == 3
        id_method = 'OE';
    end
    for seed = seed_g
        n = network_swing_simple(4, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
        n.Adj_ref = n.Adj_ref*0;
        final_model_set = cell(ITR,1);
        init_sys_set = cell(ITR,1);
        parfor itr = 1 : ITR
            try
                [final_model,init_sys] = SPEM_test_code( n, c_n, n_n, id_in_p, noise_p, id_method, sim_seed(itr));
                final_model_set(itr) = {final_model};
                init_sys_set{itr} = {init_sys};
            catch
                final_model_set(itr) = {[]};
                init_sys_set{itr} = {[]};
            end
        end
        [~, sys_env] = n.get_sys_local(c_n);
        save(strcat(location,'\',id_method,'_init_Seed_',num2str(seed)),'sys_env','final_model_set','init_sys_set');
        clear final_model_set init_sys_set sys_env
    end
end
