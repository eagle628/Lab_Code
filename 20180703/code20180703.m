% Environment Identification
%{
Localには，同定用入力が入ると仮定して，Identificationを行うとする．
Environmentに，雑音が入ったときにどのような結果を示すかを確認する．
%}
%% initiallize workspace
clear 
close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 4;
seed = 10;% Line:2, Tri:4, Square:8, Star:10
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;
n.plot()
%% Set I/O node 
c_n = 1;
n_n = [4];
%% Signal Power
id_in_p = 1;
noise_p = [0.1,0.1,0.1];
%% Identification Method
id_method = 'ARMAX';
ITR = 10;
%% Figure Name
ob_v  = {strcat('v_node',num2str(c_n))};
ob_w  = {strcat('w_node',num2str(c_n))};
ob_y  = {strcat('y_node',num2str(c_n))};
ID_in = {strcat('d_node',num2str(c_n))};
Noise = cell(1,numel(n_n));
for i = 1 : numel(n_n)
    Noise(i) = {strcat('d_node',num2str(n_n(i)))};
end

name_c = '';
for i = 1 : numel(c_n)
    name_c = strcat(name_c,num2str(c_n(i)),'+');
end
name_n = '';
for i = 1 : numel(n_n)
    name_n = strcat(name_n,num2str(n_n(i)),'+');
end
str = replace(num2str(noise_p),'.','%');
name = strcat(id_method,'_N',num2str(Node_number),'_Cnode',name_c,'_Nnode',name_n,'_IDinput',num2str(id_in_p),'_Noise',str);
%% add I/O port for identificaiton
sys_org = n.get_sys();
for i =  1 : n.N
    sys_org = n.add_io(sys_org,i, strcat('node',num2str(i)));
end

min = -5;
max =  5;
[sys_local, sys_env] = n.get_sys_local(c_n);
[mag_env,phase_env,omega] = bode(sys_env,{10^min,10^max});
mag_env=reshape(mag_env(1,1,:),[],1);
phase_env=reshape(phase_env(1,1,:),[],1);

%% 
rng(28);
rand_s = randi(1000,2,ITR);
Mag = zeros(length(omega),ITR);
Pha = zeros(length(omega),ITR);

% Model dimension
dim = 6;
% OE option
opt = oeOptions;
opt.Display = 'on';
opt.SearchOption.MaxIter = 30;
opt.SearchOption.Tolerance = 0.01;
opt.SearchMethod = 'lm';
%
%n.controllers = {};
%n.add_controller(n_n(1), diag([1,1000]), 1)
%sys_org2 = n.get_sys_controlled(sys_org);
sys_org2 = sys_org;

model_set = cell(1,ITR);
for i = 1:ITR
    % ID input
    N = 100000;
    Ts = 0.01;
    t = (0:N-1)'*Ts;
    cn1 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(c_n),'RandomStream','mt19937ar with seed','Seed',rand_s(1,i));
    d1 = cn1();
    %d1(:,1) = zeros(N,1);
    d1 = d1*id_in_p;
    cn2 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2*numel(n_n),'RandomStream','mt19937ar with seed','Seed',rand_s(2,i));
    d2 = cn2();
    for k = 1 : numel(n_n)
        %d2(:,2*k-1) = zeros(N,1);
        d2 = d2*noise_p(k);
    end
    % simulation
    lsim_type = 'foh';
    % Response of v&w 
    v = lsim(sys_org2(ob_v, cat(2,ID_in,Noise)), [d1,d2], t,lsim_type);
    w = lsim(sys_org2(ob_w, cat(2,ID_in,Noise)), [d1,d2], t,lsim_type);
    data_o = iddata(v,w,Ts);
    if string(id_method) == 'OE'
        R = 10;
        data = resample(data_o,1,R);
        init_sys = oe(data,[dim,dim,0],opt);
    elseif string(id_method) == 'ARMAX'
        R = 10;
        data = resample(data_o,1,R);
        %init_sys = armax(data,[dim,dim+1,dim,0],opt);
        init_sys = armax(data,[dim,dim+1,dim,0],opt);
    elseif string(id_method) == 'BJ'
        R = 1;
        data = resample(data_o,1,R);
        init_sys = bj(data,[dim+1,dim,dim,dim,0],opt);
    end
    model_p = d2c(init_sys);
    model = ss(model_p);
    model_set(i) = {model};
    
    [mag,phase] = bode(model,omega);
    Mag(:,i) =reshape(mag(1,1,:),[],1);
    Pha(:,i) =reshape(phase(1,1,:),[],1);
end

%% Bode Plot
% Bode Prepare
fig1 = figure('Name',strcat('Bode_',name));
ax1 = subplot(2,1,1);
hold on;
ax2 = subplot(2,1,2);
hold on;

ax1.XScale = 'log';
ax1.XMinorGrid = 'on';
ax1.YMinorGrid = 'on';
ax1.Box = 'on';
ax1.XLabel.String = 'Frequency [rad/s]';
ax1.YLabel.String = 'Gain [dB]';

ax2.XScale = 'log';
ax2.XMinorGrid = 'on';
ax2.YMinorGrid = 'on';
ax2.Box = 'on';
ax2.XLabel.String = 'Frequency [rad/s]';
ax2.YLabel.String = 'Phase [deg]';

%%%%%%%%%%%%%%
plot(ax1,omega,mag2db(mag_env),'r','linewidth',3);
plot(ax2,omega,phase_env,'r','linewidth',3);
for i = 1 : ITR
    plot(ax1,omega,mag2db(Mag(:,i)),'b');
    plot(ax2,omega,Pha(:,i),'b');
end
%% Add Controller (Extend Retro fit)
Q = kron(eye(1*numel(c_n)),diag([1,1000]));
R = kron(eye(1*numel(c_n)),diag([1e-3]));
% Ideal 
n.controllers = {};
n.add_controller(c_n, Q, R)
%n.add_controller(c_n, Q, R);
controlled_sys_I = n.get_sys_controlled(sys_org);
controlled_sys_I_K = n.controllers{1}.sys_K;

% Identification
controlled_sys_set = cell(1,ITR);
controlled_sys_K_set = cell(1,ITR); 
for i = 1 : ITR
    n.controllers = {};
    n.add_controller(c_n,model_set{i}, Q, R);
    controlled_sys_set(i) = {n.get_sys_controlled(sys_org)};
    controlled_sys_K_set(i) = {n.controllers{1}.sys_K};
end

%% Response simulation
s_t = 600;%simulation_time
time = 0:0.01:s_t-0.01;
time = time';
% Original Response
yO = impulse(sys_org(ob_y,ID_in), time);
%yO = impulse(sys_org(ob_y,Noise), time);
% Ideal Response
yI = impulse(controlled_sys_I(ob_y,ID_in), time);
%yI = impulse(controlled_sys_I(ob_y,Noise), time);
yI_de = impulse(controlled_sys_I_K('y','x'),time);

y_id_set = cell(1,ITR);
y_id_de_set = cell(1,ITR); 
for i = 1 : ITR
    sys_re = controlled_sys_set{i};
    sys_de = controlled_sys_K_set{i};
    y_id_set(i) = {impulse(sys_re(ob_y,ID_in), time)};
    %y_id_set(i) = {impulse(sys_re(ob_y,Noise), time)};
    y_id_de_set(i) = {impulse(sys_de('y','x'),time)};
end

%% Response plot
check = 1;
state = 2;
fig2 = figure('Name',strcat('Response_',name));
ax1 = subplot(3,1,1);
hold on;
ori = plot(time,yO(:,state,check),'g','linewidth',1.0);
ideal = plot(time,yI(:,state,check),'r','linewidth',1.5);
for i = 1 : ITR
    id = plot(time,y_id_set{i}(:,state,check),'b');
end
legend([ori,ideal,id],'original','Ideal','Identification');
ax2 = subplot(3,1,2);
hold on
plot(time,yI_de(:,state,check),'r','linewidth',1.5);
for i = 1 : ITR
    plot(time,y_id_de_set{i}(:,state,check),'b');
end
ax3 = subplot(3,1,3);
hold on
plot(time,yI(:,state,check)-yI_de(:,state,check),'r','linewidth',1.5);
for i = 1 : ITR
    plot(time,y_id_set{i}(:,state,check)-y_id_de_set{i}(:,state,check),'b');
end

ax1.Title.Interpreter = 'latex';
ax1.Title.FontSize = 12;
ax1.Title.String = '$y$';
ax1.Box = 'on';
ax2.Title.Interpreter = 'latex';
ax2.Title.FontSize = 12;
ax2.Title.String = '$\hat{y}$';
ax2.Box = 'on';
ax3.Title.FontSize = 12;
ax3.Title.Interpreter = 'latex';
ax3.Title.String = '$\tilde{y}$';
ax3.Box = 'on';

%% save
%savefig(fig1,fig1.Name,'compact')
%savefig(fig2,fig2.Name,'compact')
%save(name)
beep
