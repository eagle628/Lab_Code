close all
clear


Node_number_g = {4,6,10,20,30};
number_add_g = {1,2};
add_initial_g = {[1,0],[0,0],[0,1]};
power_g = {0,1};
simu_noise_g ={0,1};
%{
Node_number_g = {10};
number_add_g = {1};
add_initial_g = {[1,0]};
power_g = {0};
simu_noise_g ={0};
%}
location = 'C:\Users\Naoya Inoue\Desktop\Test\Add_node';

i = 1;
j = 1;
k = 1;
l = 1;
m = 1;

for i = 1 : numel(Node_number_g)
    for j = 1 : numel(number_add_g)
        for k = 1 : numel(add_initial_g)
            for l = 1 : numel(power_g)
                for m = 1 : numel(simu_noise_g)
                    disp('Start___________________________________')
                    name = code20180711(Node_number_g{i},number_add_g{j},add_initial_g{k},power_g{l},simu_noise_g{m},location);
                    disp(strcat('End____',name,'_______________'))
                end
            end
        end
    end
end
