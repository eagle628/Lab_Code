% Environment Identification
% other point in noise

% まずは，使えるidentification Methodを探すこと
%% initiallize workspace
clear 
close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 30;
seed = 10;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;

%% Set I/O node 
node = [1,2,3,4];

%% add I/O port for identificaiton
sys_org = n.get_sys();
for i =  1 : numel(node)
    sys_org = n.add_io(sys_org,node(i), strcat('Id_node',num2str(node(i))));
end

%% Set Identificaiton Input
N = 50000;
Ts = 0.25;
%d = randn(N, 10);
t = (0:N-1)'*Ts;
noise_power = 1;
cn = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',10);
d= cn()*noise_power;

%% Response simulation
[~, sys_env] = n.get_sys_local(1);
sys_initial = balred(sys_env, 6);

lsim_type = 'foh';

% Response of v&w for d1_node
v0 = lsim(sys_org({'v_Id_node1'}, {'d_Id_node1'}), d(:,1:4), t,lsim_type);
w0 = lsim(sys_org({'w_Id_node1'}, {'d_Id_node1'}), d(:,1:4), t,lsim_type);
    %model0 = identification(t,w0,v0,ss(model_armax0));
    %model0 = identification(Ts,w0,v0);
dim = 6;
na = dim;nb = dim+1;nc = dim;nd = dim;nf = dim;nk = 0;
nn = [na,nb,nc,nd,nf,nk];
data0 = iddata(v0,w0,Ts);
model1 = pem(data0,nn,'Focus','Simulation');
    %model0 = n4sid(data0,6,'Feedthrough',1);
model1_ss = minreal(ss(model1));


%% Chracter plot Area
figure
bode(sys_env);
hold on; grid on;
bode(sys_initial);
bode(model1);


%% Add Controller (Extend Retro fit)
c_n = 1;
sys_all = n.get_sys();
sys_all = n.add_io(sys_all,c_n, strcat('node',num2str(c_n)));
Q = diag([1,1000]);
R = diag([1e-3]);
% Obj.controllers initialize
% Retro
n.controllers = {};
n.add_controller(c_n, Q, R);
sys_cr = n.get_sys_controlled(sys_all);
% Ideal Extend retro
n.controllers = {};
n.add_controller(c_n,sys_initial, Q, R);
sys_cI = n.get_sys_controlled(sys_all);
% w/o noise
n.controllers = {};
n.add_controller(c_n,model1_ss, Q, R);
sys_c1 = n.get_sys_controlled(sys_all);

%% response
s_t = 300;%simulation_time
time = 0:0.01:s_t-0.01;
time = time';
%y = impulse(sys_all({'y_node1'}, {'d_node1'}), time);
%yr = impulse(sys_cr({'y_node1'}, {'d_node1'}), time);
%y0 = impulse(sys_c0({'y_node1'}, {'d_node1'}), time);
%y1 = impulse(sys_c1({'y_node1'}, {'d_node1'}), time);

noise = zeros(length(time),2);
noise(:,1) = randn(length(time),1);
y = lsim(sys_all({'y_node1'}, {'d_node1'}), noise, time);
yr = lsim(sys_cr({'y_node1'}, {'d_node1'}), noise, time);
yI = lsim(sys_cI({'y_node1'}, {'d_node1'}), noise, time);
y1 = lsim(sys_c1({'y_node1'}, {'d_node1'}), noise, time);


a = 4;
check = 1;
norm_number = 2;
figure
subplot(a,1,1)
plot(time,y(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y(:,1,check),norm_number)),num2str(norm(y(:,2,check),norm_number)),'location','best')
subplot(a,1,2)
plot(time,yr(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yr(:,1,check),norm_number)),num2str(norm(yr(:,2,check),norm_number)),'location','best')
subplot(a,1,3)
plot(time,yI(:,:,check),'LineWidth',1.5)
legend(num2str(norm(yI(:,1,check),norm_number)),num2str(norm(yI(:,2,check),norm_number)),'location','best')
subplot(a,1,4)
plot(time,y1(:,:,check),'LineWidth',1.5)
legend(num2str(norm(y1(:,1,check),norm_number)),num2str(norm(y1(:,2,check),norm_number)),'location','best')

%% function
function model = identification_k(t,w,v,initial)
    n_model = 6;
    m = model_cont_std(n_model, 0);%何次モデルにするか
    m.set_tf(initial);% 低次元化パラメータを初期値にセット
    m.lsim_type = 'foh';
    m.fit_mymqt(t, w, v);% マルカール法によるパラメータ推定
    model = minreal(ss(m));% 最小実現
end

function model = identification(Ts,w,v)
    % Set option
    opt = polyestOptions;
    opt.Display = 'on';
    opt.SearchMethod = 'lm';
    %opt.SearchOptions.MaxIterations = 1000;
    opt.Focus = 'prediction';
    % dcimation
    r = 1;
    v = decimate(v,r);
    w = decimate(w,r);
    data = iddata(v,w,Ts*r);
    model_dim = 12;
    na = model_dim;nb = model_dim+1;nc = model_dim;nd = model_dim;nf = model_dim;nk = 0;
        %model_ini = arx(data,[na nb nk],opt);
    model_ini = armax(data,[na nb nc nk],opt);
        %model_ini = bj(data,[nb nc nd nf nk],'Display','on','Focus','simulation');
    model = pem(data,model_ini,opt);
    model = minreal(ss(model));
    %model = balred(model,6);
end