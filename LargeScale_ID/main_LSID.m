close all
clear 
% %% node3
% root = {
%         'C:\Users\NaoyaInoue\Desktop\figure_set\node3_confirm_ss_oe_spem\node1_100\np0',...
%         'C:\Users\NaoyaInoue\Desktop\figure_set\node3_confirm_ss_oe_spem\node1_100\np001'...
%         'C:\Users\NaoyaInoue\Desktop\figure_set\node3_confirm_ss_oe_spem\node1_100\np01',...
%         };
% 
% max_iter = 100;
% n_n = [2,3];
% Node_number = 3;
% noise_p_g = [0,0.01,0.1];
% controller_number = [];
% c_n = 1;
% 
% idx = 1;
% try
%     for noise_p = noise_p_g
%         LSID_func2(Node_number, noise_p, n_n, root{idx}, max_iter, controller_number, c_n);
%         idx = idx +1 ;
%     end
%     mail_message(strcat('node3_np',num2str(idx-1)))
% catch
%     mail_message('Error')
% end

%% node30
root = {
        'C:\Users\NaoyaInoue\Desktop\figure_set\node30_confirm_ss_oe_spem\node1\np0',...
        'C:\Users\NaoyaInoue\Desktop\figure_set\node30_confirm_ss_oe_spem\node1\np001'...
        'C:\Users\NaoyaInoue\Desktop\figure_set\node30_confirm_ss_oe_spem\node1\np01',...
        };

max_iter = 100;
n_n = [2:30];
Node_number = 30;
noise_p_g = [0,0.01];
controller_number = [];
c_n = 1;

idx = 1;
try
    for noise_p = noise_p_g
        LSID_func2(Node_number, noise_p, n_n, root{idx}, max_iter, controller_number, c_n);
        idx = idx +1 ;
    end
    mail_message(strcat('node30_np',num2str(idx-1)))
catch
    mail_message('Error')
end
