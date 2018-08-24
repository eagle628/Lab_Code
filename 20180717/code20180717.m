% Environment Identification
%{
Nodeの追加したときの雑音に対する影響
retro制御器を新たに増設したとき
%}
%% initiallize workspace
%clear 
%close all
%% Function
%function name = code20180717(Node_number,control_node_a,noise_point,noise_power,location)
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 4;
seed = 10;
m = [1,2];
d = [2,10]*1e-2;
b = 1;
Y = [1,5];
%Y = [0,1];
r_c = 0.1;
n_org = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_org.Adj_ref = n_org.Adj_ref*0;

%% Controlled node number
control_node_o = 1;
%% Identification SN
id_SN = 0.1;
id_noise_node = [4];      % identification noise point
%% Identification Method
id_method = 'ARMAX'; % ARMAX or OE
%% Add node Config
number_add = 1; % number of add node
add_initial = [1,0];
%% Add Controller point
control_node_a = 2;
%% Number of Iteration
ITR = 5; 
%% simulation noise point & power
noise_point = 4;
noise_power = 0.1;
%% Save Directory
location = 'C:\Users\Naoya Inoue\Desktop\Test\add_node&controller\noise_E';
%% Node & Figure Name
% Named nodes
ob_v_p  = {strcat('v_node',num2str(control_node_o))};
ob_w_p  = {strcat('w_node',num2str(control_node_o))};
ob_y_p  = {strcat('y_node',num2str(control_node_o))};
ID_in_p = {strcat('d_node',num2str(control_node_o))};
ob_xhat_p = {strcat('xhat_controlled',num2str(control_node_o))};
con_x_p = {strcat('controller_x',num2str(control_node_o))};
id_noise_p = cell(1,numel(id_noise_node));
for i = 1 : numel(id_noise_node)
    id_noise_p(i) = {strcat('d_node',num2str(id_noise_node(i)))};
end
sim_noise_p = cell(1,numel(noise_point));
for i = 1 :  numel(noise_point)
    sim_noise_p(i) ={strcat('d_node',num2str(noise_point(i)))};
end
% Named experiment condition
name_c = '';
for i = 1 : numel(control_node_o)
    name_c = strcat(name_c,num2str(control_node_o(i)),'+');
end
name_n = '';
for i = 1 : numel(id_noise_node)
    name_n = strcat(name_n,num2str(id_noise_node(i)),'+');
end
str1 = replace(num2str(id_SN),'.','%');
str2 = replace(num2str(noise_power),'.','-');
%
name = strcat('N',num2str(Node_number),'_Cnode',num2str(control_node_o),'_ACnode',num2str(control_node_a),...
                '_NoisePower',str2,'_NoiseLocation_',num2str(noise_point)...
                );
%}
%% add I/O port for original system
sys_org = n_org.get_sys();
for i =  1 : numel(control_node_o)
    sys_org = n_org.add_io(sys_org,control_node_o(i), strcat('node',num2str(control_node_o(i))));
end
for i =  1 : numel(id_noise_node)
    sys_org = n_org.add_io(sys_org,id_noise_node(i), strcat('node',num2str(id_noise_node(i))));
end
for i =  1 : numel(noise_point)
    sys_org = n_org.add_io(sys_org,id_noise_node(i), strcat('node',num2str(noise_point(i))));
end
%% Generate Add node System 
% Generate network for Add node
n_add = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_add.Adj_ref = n_add.Adj_ref*0;
% Make New node
for i = 1 : number_add
    n_add.add_node(m, d, b, Y, r_c);
end
% New Network add I/O port
sys_add = n_add.get_sys();
for i =  1 : numel(control_node_o)
    sys_add = n_add.add_io(sys_add,control_node_o(i), strcat('node',num2str(control_node_o(i))));
end
for i =  1 : numel(id_noise_node)
    sys_add = n_add.add_io(sys_add,id_noise_node(i), strcat('node',num2str(id_noise_node(i))));
end
for i =  1 : numel(noise_point)
    sys_add = n_add.add_io(sys_add,id_noise_node(i), strcat('node',num2str(noise_point(i))));
end
%% Environment Character (Original)
% Chose Range Fruquency
min = -6;
max =  6;
[~, sys_env_h] = n_org.get_sys_local(control_node_o);
sys_env_o1 = balred(sys_env_h,6);
[mag_env,pha_env,omega] = bode(sys_env_h,{10^min,10^max});
mag_env=squeeze(mag_env(1,1,:));
pha_env=squeeze(pha_env(1,1,:));
% Environment for add controller
[~, sys_env_o2] = n_org.get_sys_local(control_node_a);
sys_env_o2 = balred(sys_env_o2,6);
%% Environment Character (add_system)
[~, sys_env_a_h] = n_add.get_sys_local(control_node_o);
sys_env_a = balred(sys_env_a_h,6);
[mag_env_a,pha_env_a,~] = bode(sys_env_a_h,omega);
mag_env_a = squeeze(mag_env_a(1,1,:));
pha_env_a = squeeze(pha_env_a(1,1,:));
%% LQR Controller Parameter
Q = kron(eye(1*numel(control_node_o)),diag([1,1]));
R = kron(eye(1*numel(control_node_o)),diag([1e-3]));
%% Add Controller (Smiple) for original system
n_org.controllers = {};
n_org.add_controller(control_node_o,ss([],[],[],0), Q, R);
%n_org.add_controller(control_node_o,sys_env_o1, Q, R);
controlled_sys_S = n_org.get_sys_controlled(sys_org);
%% Add 2 Retro Controller for original system (extend & simple)
n_org.controllers = {};
n_org.add_controller(control_node_o,ss([],[],[],0), Q, R);
%n_org.add_controller(control_node_o,sys_env_o1, Q, R);
n_org.add_controller(control_node_a,ss([],[],[],0), Q, R);
%n_org.add_controller(control_node_a,sys_env_o2, Q, R);
controlled_sys_a_c_S = n_org.get_sys_controlled(sys_org);
%% Add Controller (Extend) for add node system
n_add.controllers = {};
n_add.add_controller(control_node_o,ss([],[],[],0), Q, R);
%n_add.add_controller(control_node_o,sys_env_o1, Q, R);
controlled_sys_a_n_S = n_add.get_sys_controlled(sys_add);
%% Identification
% Initial And Get memory
rng(6);
rand_s = randi(1000,2,ITR);
Mag = zeros(length(omega),ITR);
Pha = zeros(length(omega),ITR);
model_set = cell(1,ITR);
% Model dimension
dim = 6;
% OE option
opt = oeOptions;
opt.Display = 'on';
opt.SearchOption.MaxIter = 100;
opt.SearchOption.Tolerance = 1e-3;
opt.SearchMethod = 'lm';
% Identification
for i = 1:ITR
    % ID input
    N = 100000;
    Ts = 0.01;
    t = (0:N-1)'*Ts;
    cn1 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(control_node_o),'RandomStream','mt19937ar with seed','Seed',rand_s(1,i));
    id_noise1 = cn1();
    clear cn1
    cn2 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(id_noise_node),'RandomStream','mt19937ar with seed','Seed',rand_s(2,i));
    id_noise2 = cn2();
    clear cn2
    id_noise2 = id_noise2*id_SN;
    % simulation
    lsim_type = 'foh';
    % Response of v&w 
    v = lsim(sys_org(ob_v_p, cat(2,ID_in_p,id_noise_p)), [id_noise1,id_noise2], t,lsim_type);
    w = lsim(sys_org(ob_w_p, cat(2,ID_in_p,id_noise_p)), [id_noise1,id_noise2], t,lsim_type);
    data_o = iddata(v,w,Ts);
    if string(id_method) == 'OE'
        data = resample(data_o,1,10);
        init_sys = oe(data,[dim,dim,0],opt);
    else
        data = resample(data_o,1,1);
        init_sys = armax(data,[dim,dim+1,dim,0],opt);
    end
    model_p = d2c(init_sys);
    model = ss(model_p);
    model_set(i) = {model};
    
    [mag,phase] = bode(model,omega);
    Mag(:,i) = squeeze(mag(1,1,:));
    Pha(:,i) = squeeze(phase(1,1,:));
end
%% Identificaiton Character
fig1 = figure_config.set_figure_bode(name);
figure_config.plot_bode(fig1.axes,omega,mag_env,pha_env,{'r-','linewidth',3});
figure_config.plot_bode(fig1.axes,omega,mag_env_a,pha_env_a,{'r:','linewidth',3});
figure_config.plot_bode(fig1.axes,omega,Mag,Pha,{'b-','linewidth',0.8});

add_system = fig1.axes.ax1.Children(end-1);
ori_system = fig1.axes.ax1.Children(end);
extend = fig1.axes.ax1.Children(1);
legend(fig1.axes.ax1,[add_system,ori_system,extend],'add','original','identificaiton','location','best')
%% Add Controller (Extend) for original system (Identification)
controlled_sys_id_set = cell(ITR,1);
for i = 1: ITR
    n_org.controllers = {};
    n_org.add_controller(control_node_o,model_set{i}, Q, R);
    controlled_sys_id_set(i) = {n_org.get_sys_controlled(sys_org)};
end
%% Add 2 Retro Controller for original system (Identification Extend & Simple)
controlled_sys_a_c_id_set = cell(ITR,1);
n.org.controllers = {};
for i = 1: ITR
    n_org.controllers = {};
    n_org.add_controller(control_node_o,model_set{i}, Q, R);        % Extend
    n_org.add_controller(control_node_a,ss([],[],[],0), Q, R);      % Simple
    %n_org.add_controller(control_node_a,sys_env_o2, Q, R);           % Extend
    controlled_sys_a_c_id_set(i) = {n_org.get_sys_controlled(sys_org)};
end
%% Add Controller for add node system (Identification Extend)
controlled_sys_a_n_id_set = cell(ITR,1);
for i = 1: ITR
    n_add.controllers = {};
    n_add.add_controller(control_node_o,model_set{i}, Q, R);
    controlled_sys_a_n_id_set(i) = {n_add.get_sys_controlled(sys_add)};
end
%% simulation
simulation_time = 100;%simulation_time
Ts_s = 0.01;
t_s = 0:Ts_s:simulation_time-Ts_s;
t_s = t_s';

sim_noise = zeros(length(t_s),4);
sim_seed = 1024;
cn3 = dsp.ColoredNoise('Color','white','SamplesPerFrame',length(t_s),'NumChannels',4,'RandomStream','mt19937ar with seed','Seed',sim_seed);
noise = cn3()*noise_power;
clear cn3
sim_noise(:,3) = noise(:,3);
%sim_noise(:,2) = noise(:,2);
sim_noise(round(length(t_s)/2),2) = 1/Ts_s;
%sim_noise(round(length(t_s)/4),1) = 1/Ts_s;

%%%%%% Original Response
y_o = lsim(sys_org( ob_y_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
y_a_n = lsim(sys_add( ob_y_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);

%%%%%% Simple retro Response (Controlled)
% Original
y_S_o = lsim(controlled_sys_S( ob_y_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
xhat_S_o = lsim(controlled_sys_S( ob_xhat_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
% Add node
y_S_a_n = lsim(controlled_sys_a_n_S( ob_y_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
xhat_S_a_n = lsim(controlled_sys_a_n_S( ob_xhat_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
% Add controller
y_S_a_c = lsim(controlled_sys_a_c_S( ob_y_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
xhat_S_a_c = lsim(controlled_sys_a_c_S( ob_xhat_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);

%%%%%% Extend Response 
% Get Memory
y_id_o = zeros(length(t_s),2,ITR);
xhat_id_o = zeros(length(t_s),2,ITR);
y_id_a_n = zeros(length(t_s),2,ITR);
xhat_id_a_n = zeros(length(t_s),2,ITR);
y_id_a_c = zeros(length(t_s),2,ITR);
xhat_id_a_c = zeros(length(t_s),2,ITR);
for i = 1 : ITR 
    % original system
    y_id_o(:,:,i) = lsim(controlled_sys_id_set{i}( ob_y_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
    xhat_id_o(:,:,i) = lsim(controlled_sys_id_set{i}( ob_xhat_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
    % add node system
    y_id_a_n(:,:,i) = lsim(controlled_sys_a_n_id_set{i}( ob_y_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
    xhat_id_a_n(:,:,i) = lsim(controlled_sys_a_n_id_set{i}( ob_xhat_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
    % add controller system
    y_id_a_c(:,:,i) = lsim(controlled_sys_a_c_id_set{i}( ob_y_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
    xhat_id_a_c(:,:,i) = lsim(controlled_sys_a_c_id_set{i}( ob_xhat_p, { ID_in_p, sim_noise_p}), sim_noise, t_s);
end

%%%%%%%% Drawing
state = 2;
fig3 = figure_config.set_figure_retro(strcat('Response_original','(',name,')'));
fig4 = figure_config.set_figure_retro(strcat('Response_add node','(',name,')'));
fig5 = figure_config.set_figure_retro(strcat('Response_add controller','(',name,')'));

% Original
figure_config.plot_retro2( fig3.axes, t_s,y_o(:,state), [], {'g-','linewidth',2.0});
figure_config.plot_retro2( fig4.axes, t_s,y_a_n(:,state), [], {'g-','linewidth',2.0});
figure_config.plot_retro2( fig5.axes, t_s,y_o(:,state), [], {'g-','linewidth',2.0});
% Simple
figure_config.plot_retro2( fig3.axes, t_s,y_S_o(:,state), xhat_S_o(:,state), {'r-','linewidth',2.0});
figure_config.plot_retro2( fig4.axes, t_s,y_S_a_n(:,state), xhat_S_a_n(:,state), {'r-','linewidth',2.0});
figure_config.plot_retro2( fig5.axes, t_s,y_S_a_c(:,state), xhat_S_a_c(:,state), {'r-','linewidth',2.0});
% Extend
figure_config.plot_retro2( fig3.axes, t_s,squeeze(y_id_o(:,state,:)), squeeze(xhat_id_o(:,state,:)), {'b-','linewidth',0.8});
figure_config.plot_retro2( fig4.axes, t_s,squeeze(y_id_a_n(:,state,:)), squeeze(xhat_id_a_n(:,state,:)), {'b-','linewidth',0.8});
figure_config.plot_retro2( fig5.axes, t_s,squeeze(y_id_a_c(:,state,:)), squeeze(xhat_id_a_c(:,state,:)), {'b-','linewidth',0.8});

% Range Limit
fig3.axes.ax2.YLim = [ -0.1, 0.1];
fig4.axes.ax2.YLim = [ -0.1, 0.1];
fig5.axes.ax2.YLim = [ -0.1, 0.1];

% Legend
ori = fig3.axes.ax1.Children(end);
simple = fig3.axes.ax1.Children(end-1);
extend = fig3.axes.ax1.Children(end-2);
legend(fig3.axes.ax1,[ori,simple,extend],'original','simple','extend')

ori = fig4.axes.ax1.Children(end);
simple = fig4.axes.ax1.Children(end-1);
extend = fig4.axes.ax1.Children(end-2);
legend(fig4.axes.ax1,[ori,simple,extend],'original','simple','extend')

ori = fig5.axes.ax1.Children(end);
simple = fig5.axes.ax1.Children(end-1);
extend = fig5.axes.ax1.Children(end-2);
legend(fig5.axes.ax1,[ori,simple,extend],'original','simple-simple','extend-simple')
%% save
%{
savefig(fig1.f,strcat(location,'\',fig1.f.Name));
savefig(fig3.f,strcat(location,'\',fig3.f.Name));
savefig(fig4.f,strcat(location,'\',fig4.f.Name));
savefig(fig5.f,strcat(location,'\',fig5.f.Name));
clear fig1 fig3 fig4 fig5 
save(name);
%}

%end