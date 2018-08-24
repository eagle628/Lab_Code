clear
close all
%% Location
location = 'C:\Users\Naoya Inoue\Desktop\Test\2018SummerReport\SPEM';
%% File Name
name = 'ARX_init_Seed_8';
%% Load Data
load(strcat(location,'\',name))
%% Environment Character (Original)
% Chose Range Fruquency
min = -4;
max =  2;
[mag_env,pha_env,omega] = bode(sys_env,{10^min,10^max});
mag_env=squeeze(mag_env(1,1,:));
pha_env=squeeze(pha_env(1,1,:));
%% size
ITR = numel(final_model_set);
%% Identification Character
Mag_Id = zeros(length(omega),ITR);
Pha_Id = zeros(length(omega),ITR);
for i = 1 : ITR
    [mag,pha] = bode(final_model_set{i},omega);
    mag = squeeze(mag(1,1,:));
    pha = squeeze(pha(1,1,:));
    Mag_Id(:,i) = mag;
    Pha_Id(:,i) = pha;
end
%% init_sys Character
Mag_init = zeros(length(omega),ITR);
Pha_init = zeros(length(omega),ITR);
for i = 1 : ITR
    [mag,pha] = bode(final_model_set{i},omega);
    mag = squeeze(mag(1,1,:));
    pha = squeeze(pha(1,1,:));
    Mag_init(:,i) = mag;
    Pha_init(:,i) = pha;
end
%% Drawing
fig1 = figure_config.set_figure_bode(name,[-600,0,300,400]);
figure_config.plot_bode(fig1.axes,omega,mag_env,pha_env,{'r-','linewidth',3});
figure_config.plot_bode(fig1.axes,omega,Mag_init,Pha_init,{'g-','linewidth',1.5});
figure_config.plot_bode(fig1.axes,omega,Mag_Id,Pha_Id,{'b-','linewidth',0.8});

Env = fig1.axes.ax1.Children(end);
Id = fig1.axes.ax1.Children(1);
init = fig1.axes.ax1.Children(end-1);
legend(fig1.axes.ax1 ,[ Env, init, Id] ,'Environment','init system','Identification','location','best')