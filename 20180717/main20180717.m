%% initiallize workspace
clear 
close all
%%
Node_number_g = [4];
%Node_number_g = [4,10,20,40];
%%
control_node_a_g = [2,4];
%%
noise_point_g = [2,4];
%%
noise_power_g = [0.1,1];
%% Save directory
location = 'C:\Users\Naoya Inoue\Desktop\Test\add_node&controller\noise_E';
%% 
ITR = 100;
%% 
for i = 1 : numel(Node_number_g)
    Node_number = Node_number_g(i);
    for j = 1 : numel(control_node_a_g)
        control_node_a = control_node_a_g(j);
        for k = 1 : numel(noise_point_g)
            noise_point = noise_point_g(k);
            for l = 1 : numel(noise_power_g)
                noise_power = noise_power_g(l);
                disp('Start______')
                name = code20180717(Node_number,control_node_a,noise_point,noise_power,location);
                disp(strcat('END______',name))
                close all
            end
        end
    end
end

beep