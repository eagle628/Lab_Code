clear
close all

load('test_figure_data')

check = 1;
state = 2;

H = figure_config.set_figure_retro();

figure_config.plot_retro(H.axes,time...
                    ,yO(:,state,check)...
                    ,[]...
                    ,{'g-','Linewidth',1.0});

figure_config.plot_retro(H.axes,time...
                    ,yI(:,state,check)...
                    ,yI_de(:,state,check)...
                    ,{'r-','Linewidth',1.5});
figure_config.plot_retro(H.axes,time...
                    ,squeeze(y_id_set(:,state,check,:))...
                    ,squeeze(y_id_de_set(:,state,check,:))...
                   ,{'b-','Linewidth',0.8});

               %
ori = H.axes.ax1.Children(end);
ideal = H.axes.ax1.Children(end-1);
ID = H.axes.ax1.Children(end-2);

legend(H.axes.ax1,[ori,ideal,ID],'original','Ideal','identificaiton')
                   %}