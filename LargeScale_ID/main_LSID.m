close all
clear 
%% node60
root = 'C:\Users\NaoyaInoue\Desktop\figure_set\LSID\node100\nocontroller';

max_iter = 20;
n_n = [2:100];
Node_number = 100;
noise_p_g = [0, 0.01, 0.1];

for noise_p = noise_p_g
    if noise_p == 0
        name = '0';
    elseif noise_p == 0.1
        name = '01';
    elseif noise_p == 0.01
        name = '001';
    elseif noise_p == 1
        name = '1';
    end
    root_dir = strcat(root,'\np',name'); 
    LSID_func(Node_number, noise_p, n_n, root_dir, max_iter);
end

