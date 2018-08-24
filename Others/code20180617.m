% Environment Identification
% 2 indetification input
%% initiallize workspace
clear 
close all
%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
Node_number = 4;
seed = 10;
n = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n.Adj_ref = n.Adj_ref*0;

%% Set I/O node 
c_n = 1;
n_n = 2;
name = strcat('N',num2str(Node_number),'_Cnode',num2str(c_n),'_Nnode',num2str(n_n));
%% add I/O port for identificaiton
sys_org = n.get_sys();
for i =  1 : Node_number
    sys_org = n.add_io(sys_org,i, strcat('node',num2str(i)));
end

[sys_local, sys_env] = n.get_sys_local(c_n);
sys_local_vw = ss(sys_local.A,sys_local.B(:,sys_local.InputGroup.v),sys_local.C(sys_local.OutputGroup.w,:),[]);


%% figure
fig1 = figure('Name',strcat(name,'_LocalID'));
hold on
fig2 = figure('Name',strcat(name,'_EnvID'));
hold on
fig3 = figure('Name',strcat(name,'_BothlID'));
hold on

ITR = 10;
rng(6);
rand_s = randi(1000,2,ITR);
for i = 1 : ITR
    %% Set Identificaiton Input
    N = 50000;
    Ts = 0.01;
    t = (0:N-1)'*Ts;
    cn1 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2,'RandomStream','mt19937ar with seed','Seed',rand_s(1,i));
    d1 = cn1();
    d1 = d1*1;
    cn2 = dsp.ColoredNoise('Color','white','SamplesPerFrame',N,'NumChannels',2,'RandomStream','mt19937ar with seed','Seed',rand_s(2,i));
    d2 = cn2();
    d2 = d2*1;

    %% Response simulation
    lsim_type = 'foh';
    R = 10;
    % Model dimension
    dim = 6;
    opt = oeOptions;
    opt.Display = 'on';
    opt.SearchOption.MaxIter = 50;%(dim*2+1)*100;
    opt.SearchOption.Tolerance = 1e-4;
    % Response of v&w in single id input (local)
    v_l = lsim(sys_org({strcat('v_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), d1, t,lsim_type);
    w_l = lsim(sys_org({strcat('w_node',num2str(c_n))}, {strcat('d_node',num2str(c_n))}), d1, t,lsim_type);
    data_l_o = iddata(v_l,w_l,Ts);
    data_l = resample(data_l_o,1,R);
    init_sys_l = oe(data_l,[dim,dim,0],opt);
    model_l_p = d2c(init_sys_l);
    model_l = ss(model_l_p);
    % Response of v&w in single id input (environment)
    v_e = lsim(sys_org({strcat('v_node',num2str(c_n))}, {strcat('d_node',num2str(n_n))}), d2, t,lsim_type);
    w_e = lsim(sys_org({strcat('w_node',num2str(c_n))}, {strcat('d_node',num2str(n_n))}), d2, t,lsim_type);
    data_e_o = iddata(v_e,w_e,Ts);
    data_e = resample(data_e_o,1,R);
    init_sys_e = oe(data_e,[dim,dim,0],opt);
    %init_sys_e = armax(data_e,[6,6,6,0],opt);
    model_e_p = d2c(init_sys_e);
    model_e = ss(model_e_p);
    % Response of v&w in double id input (local&environment)
    v_le = lsim(sys_org({strcat('v_node',num2str(c_n))}, {strcat('d_node',num2str(c_n)),strcat('d_node',num2str(n_n))}), [d1,d2], t,lsim_type);
    w_le = lsim(sys_org({strcat('w_node',num2str(c_n))}, {strcat('d_node',num2str(c_n)),strcat('d_node',num2str(n_n))}), [d1,d2], t,lsim_type);
    data_le_o = iddata(v_le,w_le,Ts);
    data_le = resample(data_le_o,1,R);
    init_sys_le = oe(data_le,[dim,dim,0],opt);
    model_le_p = d2c(init_sys_le);
    model_le = ss(model_le_p);

    %% character area
    figure(fig1);
    bode(model_l,'r')
    figure(fig2);
    bode(model_e,'b')
    figure(fig3);
    bode(model_le,'g')

end

%% 
figure(fig1);
bode(sys_env,'k:')
figure(fig2);
bode(sys_env,'k:')
figure(fig3);
bode(sys_env,'k:')

%% 
%{
saveas(fig1,fig1.Name,'m')
saveas(fig2,fig2.Name,'m')
saveas(fig3,fig3.Name,'m')
%}