close all
clear all
ITR = 100;
y1 = zeros(60000,2,ITR);
y2 = zeros(60000,2,ITR);
y3 = zeros(60000,2,ITR);
y4 = zeros(60000,2,ITR);
xhat1 = zeros(60000,2,ITR);
xhat2 = zeros(60000,2,ITR);
xhat3 = zeros(60000,2,ITR);
xhat4 = zeros(60000,2,ITR);

for i = 1 : ITR
    [y1(:,:,i),y2(:,:,i),y3(:,:,i),y4(:,:,i),xhat1(:,:,i),xhat2(:,:,i),xhat3(:,:,i),xhat4(:,:,i),t_s] = code20180724_2();
end

y1 = mean(y1,3);
y2 = mean(y2,3);
y3 = mean(y3,3);
y4 = mean(y4,3);
xhat1 = mean(xhat1,3);
xhat2 = mean(xhat2,3);
xhat3 = mean(xhat3,3);
xhat4 = mean(xhat4,3);

fig1 = figure_config.set_figure_retro('Response');
state = 2;
% Not Changes
figure_config.plot_retro2( fig1.axes, t_s,y3(:,state), xhat3(:,state), {'g-','linewidth',3.0});
% Add Node
figure_config.plot_retro2( fig1.axes, t_s,y1(:,state), xhat1(:,state), {'b-','linewidth',0.8});
% Add Controller
figure_config.plot_retro2( fig1.axes, t_s,y2(:,state), xhat2(:,state), {'r-','linewidth',0.8});
% Remove Controller
figure_config.plot_retro2( fig1.axes, t_s,y4(:,state), xhat4(:,state), {'k-','linewidth',0.8});

ex2_ori = fig1.axes.ax1.Children(end);
ex2_add = fig1.axes.ax1.Children(end-1);
ex3_ori = fig1.axes.ax1.Children(end-2);
ex1_ori = fig1.axes.ax1.Children(end-3);
legend(fig1.axes.ax1,[ex2_ori,ex2_add,ex3_ori,ex1_ori],'Original(Extend2)','AddNode(Extend2)','Original(Extend3)','Original(Extend1)','location','best')