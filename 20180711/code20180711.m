% Environment Identification
%{
Nodeを追加したときの挙動
%}

function name = code20180711(Node_number,number_add,add_initial,power,simu_noise,location)
%% initiallize workspace
%clear 
%close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
%Node_number = 4;
seed = 10;
m = [1,2];
d = [2,10]*1e-2;
b = 1;
Y = [1,5];
r_c = 0.1;
n_org = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_org.Adj_ref = n_org.Adj_ref*0;

%% Set I/O node 
c_n = 1;
n_n = [2];
%% Signal Power
id_in_p = 1;
noise_p = 0.01;
%% Identification Method
id_method = 'ARMAX';
%% Add node Config
%number_add = 1;
%add_initial = [1,0];
%% Iteration
ITR = 20; 
%% simulation noise of original system before adding node
%power = 0;
%simu_noise = 1; % Environmet noise Or Local Noise
%% Figure Name
% Named nodes
ob_v  = {strcat('v_node',num2str(c_n))};
ob_w  = {strcat('w_node',num2str(c_n))};
ob_y  = {strcat('y_node',num2str(c_n))};
ID_in = {strcat('d_node',num2str(c_n))};
ob_xhat = {strcat('xhat_controlled',num2str(c_n))};
con_x = {strcat('controller_x',num2str(c_n))};
Noise = cell(1,numel(n_n));
for i = 1 : numel(n_n)
    Noise(i) = {strcat('d_node',num2str(n_n(i)))};
end
% Named experiment condition
name_c = '';
for i = 1 : numel(c_n)
    name_c = strcat(name_c,num2str(c_n(i)),'+');
end
name_n = '';
for i = 1 : numel(n_n)
    name_n = strcat(name_n,num2str(n_n(i)),'+');
end
str1 = replace(num2str(noise_p),'.','%');
if simu_noise ==1
    str2 = 'Local';
else
    str2 = 'Env';
end
name = strcat(id_method,'_N',num2str(Node_number),'_Cnode',name_c,'_AddN',num2str(number_add),...
                '_Addinitial_',num2str(add_initial(1)),'_',num2str(add_initial(2)),...
                '_NoisePower',num2str(power),'_NoiseLocation_',str2...
                );
%% add I/O port for identificaiton
sys_org = n_org.get_sys();
for i =  1 : numel(c_n)
    sys_org = n_org.add_io(sys_org,c_n(i), strcat('node',num2str(c_n(i))));
end
for i =  1 : numel(n_n)
    sys_org = n_org.add_io(sys_org,n_n(i), strcat('node',num2str(n_n(i))));
end
%% Add node
% Generate network for Add node
n_add = network_swing_simple(Node_number, m, d, b, Y, r_c, seed);
n_add.Adj_ref = n_add.Adj_ref*0;
% Make New node
for i = 1 : number_add
    n_add.add_node(m, d, b, Y, r_c);
end
% New Network SS
sys_add = n_add.get_sys();
for i =  1 : numel(c_n)
    sys_add = n_add.add_io(sys_add,c_n(i), strcat('node',num2str(c_n(i))));
end
for i =  1 : numel(n_n)
    sys_add = n_add.add_io(sys_add,n_n(i), strcat('node',num2str(n_n(i))));
end
%% Environment Character (Original)
min = -10;
max =  10;
[~, sys_env_h] = n_org.get_sys_local(c_n);
sys_env = balred(sys_env_h,6);
[mag_env,pha_env,omega] = bode(sys_env_h,{10^min,10^max});
mag_env=squeeze(mag_env(1,1,:));
pha_env=squeeze(pha_env(1,1,:));
%% Environment Character (add_system)
[~, sys_env_a_h] = n_add.get_sys_local(c_n);
sys_env_a = balred(sys_env_a_h,6);
[mag_env_a,pha_env_a,~] = bode(sys_env_a_h,omega);
mag_env_a = squeeze(mag_env_a(1,1,:));
pha_env_a = squeeze(pha_env_a(1,1,:));
%% Add Controller (Extend) for original system (Ideal)
Q = kron(eye(1*numel(c_n)),diag([1,1000]));
R = kron(eye(1*numel(c_n)),diag([1e-3]));
% Ideal 
n_org.controllers = {};
n_org.add_controller(c_n,sys_env, Q, R);
controlled_sys_I = n_org.get_sys_controlled(sys_org);
controlled_sys_I_K = n_org.controllers{1}.sys_K;
%% Add Controller (Extend) for add node system (Ideal)
% Ideal 
n_add.controllers = {};
n_add.add_controller(c_n,sys_env, Q, R);
controlled_sys_add_I = n_add.get_sys_controlled(sys_add);
controlled_sys_add_I_K = n_add.controllers{1}.sys_K;
%% Identification
% Initial And Get memory
rng(28);
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
    cn1 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(c_n),'RandomStream','mt19937ar with seed','Seed',rand_s(1,i));
    d1 = cn1();
    d1 = d1*id_in_p;
    cn2 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(n_n),'RandomStream','mt19937ar with seed','Seed',rand_s(2,i));
    d2 = cn2();
    d2 = d2*noise_p;
    % simulation
    lsim_type = 'foh';
    R = 1;
    % Response of v&w 
    v = lsim(sys_org(ob_v, cat(2,ID_in,Noise)), [d1,d2], t,lsim_type);
    w = lsim(sys_org(ob_w, cat(2,ID_in,Noise)), [d1,d2], t,lsim_type);
    data_o = iddata(v,w,Ts);
    data = resample(data_o,1,R);
    if string(id_method) == 'OE'
        init_sys = oe(data,[dim,dim,0],opt);
    else
        init_sys = armax(data,[dim,dim+1,dim,0],opt);
        %init_sys = myARX(w,v,[dim,dim+1,0],Ts);
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
ID = fig1.axes.ax1.Children(1);
legend(fig1.axes.ax1,[add_system,ori_system,ID],'add','original','identificaiton','location','best')
%% Add Controller (Extend) for original system (Identification)
controlled_sys_id_set = cell(ITR,1);
controlled_sys_id_K_set = cell(ITR,1);
for i = 1: ITR
    n_org.controllers = {};
    n_org.add_controller(c_n,model_set{i}, Q, R);
    controlled_sys_id_set(i) = {n_org.get_sys_controlled(sys_org)};
    controlled_sys_id_K_set(i) = {n_org.controllers{1}.sys_K};
end
%% Add Controller (Extend) for add node system (Identification)
controlled_sys_add_id_set = cell(ITR,1);
controlled_sys_add_id_K_set = cell(ITR,1);
for i = 1: ITR
    n_add.controllers = {};
    n_add.add_controller(c_n,model_set{i}, Q, R);
    controlled_sys_add_id_set(i) = {n_add.get_sys_controlled(sys_add)};
    controlled_sys_add_id_K_set(i) = {n_add.controllers{1}.sys_K};
end
%% response simulation
% the condition of Environmet disturbance is same.
% So, Original = Ideal = Identificaiton
state = 2;
%%%%%%%%%%%%%%%%
simulation_time = 600;%simulation_time
Ts_s = 0.01;
t_s = 0:Ts_s:simulation_time-0.01;
t_s = t_s';

u = rand(length(t_s),1)*power;
u = [zeros(length(u),1),u];

add_initial = kron(ones(1,number_add),add_initial);

if simu_noise == 1
    y_o = simulation_change_system(sys_org,sys_add,u,t_s,ob_y,ID_in,add_initial);
    [y_e,v_e,w_e] = simulation_change_system(controlled_sys_I,controlled_sys_add_I,u,t_s,ob_y,ID_in,add_initial,ob_xhat,con_x);
    y_id = zeros(length(y_e),ITR,2);
    v_id = zeros(length(v_e),ITR);
    w_id = zeros(length(w_e),ITR);
    for i = 1 : ITR
        [y_id(:,i,:),v_id(:,i),w_id(:,i)] = simulation_change_system(controlled_sys_id_set{i},controlled_sys_add_id_set{i},u,t_s,ob_y,ID_in,add_initial,ob_xhat,con_x);
    end
else
    y_o = simulation_change_system(sys_org,sys_add,u,t_s,ob_y,Noise,add_initial,[],[],50);
    [y_e,v_e,w_e] = simulation_change_system(controlled_sys_I,controlled_sys_add_I,u,t_s,ob_y,Noise,add_initial,ob_xhat,con_x,50);
    y_id = zeros(length(y_e),ITR,2);
    v_id = zeros(length(v_e),ITR);
    w_id = zeros(length(w_e),ITR);
    for i = 1 : ITR
        [y_id(:,i,:),v_id(:,i),w_id(:,i)] = simulation_change_system(controlled_sys_id_set{i},controlled_sys_add_id_set{i},u,t_s,ob_y,Noise,add_initial,ob_xhat,con_x,50);
    end
end

fig2 = figure('Name',name,'position',[0,0,600,900]);
subplot(3,1,1)
plot(t_s,y_o(:,state),'g-','linewidth',0.8);
title('Original')
subplot(3,1,2)
plot(t_s,y_e(:,state),'r-','linewidth',0.8);
title('Ideal')
subplot(3,1,3)
plot(t_s,y_id(:,:,state),'b-','linewidth',0.8);
title('Identification')

%% Node add Shock identification
% The condition of Environmenet Impluse Noise is same.

% Original Response
y_o = impulse(sys_org(ob_y,ID_in), t_s);
y_a = impulse(sys_add(ob_y,ID_in), t_s);

% Ideal Response (Controlled)
y_I_o = impulse(controlled_sys_I(ob_y,ID_in), t_s);
y_I_de_o = impulse(controlled_sys_I_K('y','x'), t_s);

y_I_a = impulse(controlled_sys_add_I(ob_y,ID_in),t_s);
y_I_de_a = impulse(controlled_sys_add_I_K('y','x'),t_s);
% Identificaiton Response (Controlled)
    % get memory
y_id_o = zeros(length(t_s),2,2,ITR);
y_id_de_o = zeros(length(t_s),2,2,ITR);
y_id_a = zeros(length(t_s),2,2,ITR);
y_id_de_a = zeros(length(t_s),2,2,ITR);
% simulation
for i = 1 : ITR 
    y_id_o(:,:,:,i) = impulse(controlled_sys_id_set{i}(ob_y,ID_in),t_s);
    y_id_de_o(:,:,:,i) = impulse(controlled_sys_id_K_set{i}('y','x'), t_s);

    y_id_a(:,:,:,i) = impulse(controlled_sys_add_id_set{i}(ob_y,ID_in),t_s);
    y_id_de_a(:,:,:,i) = impulse(controlled_sys_add_id_K_set{i}('y','x'), t_s);
end
%% Response Drawing
check = 1;
state = 2;
fig4 = figure_config.set_figure_retro(strcat('Response_ori_',name));
fig3 = figure_config.set_figure_retro(strcat('Response_add_',name));

% Original
figure_config.plot_retro(fig4.axes,t_s,y_o(:,state,check),[],{'g-','linewidth',2.0});
figure_config.plot_retro(fig3.axes,t_s,y_a(:,state,check),[],{'g-','linewidth',2.0});
% Ideal
figure_config.plot_retro(fig4.axes,t_s,y_I_o(:,state,check),y_I_de_o(:,state,check),{'r-','linewidth',3});
figure_config.plot_retro(fig3.axes,t_s,y_I_a(:,state,check),y_I_de_a(:,state,check),{'r-','linewidth',3});
% Identification
figure_config.plot_retro(fig4.axes,t_s...
                        ,reshape(y_id_o(:,state,check,:),length(t_s),ITR)...
                        ,reshape(y_id_de_o(:,state,check,:),length(t_s),ITR)...
                        ,{'b-','linewidth',0.5});
figure_config.plot_retro(fig3.axes,t_s...
                        ,reshape(y_id_a(:,state,check,:),length(t_s),ITR)...
                        ,reshape(y_id_de_a(:,state,check,:),length(t_s),ITR)...
                        ,{'b-','linewidth',0.5});

%}
ori = fig4.axes.ax1.Children(end);
ideal = fig4.axes.ax1.Children(end-1);
ID = fig4.axes.ax1.Children(end-2);
legend(fig4.axes.ax1,[ori,ideal,ID],'original','Ideal','identificaiton')

ori = fig3.axes.ax1.Children(end);
ideal = fig3.axes.ax1.Children(end-1);
ID = fig3.axes.ax1.Children(end-2);
legend(fig3.axes.ax1,[ori,ideal,ID],'original','Ideal','identificaiton')


%% Save
savefig(fig1.f,strcat(location,'\',fig1.f.Name))
savefig(fig2,strcat(location,'\','Response(AddNodeShock)_',fig2.Name))
savefig(fig3.f,strcat(location,'\',fig3.f.Name))
savefig(fig4.f,strcat(location,'\',fig4.f.Name))
close all
clear fig1 fig2 fig3 fig4
save(strcat(location,'\',name));
%% For drawing Function


end