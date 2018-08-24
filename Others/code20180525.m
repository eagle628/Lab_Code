% このコードはレトロフィット制御器を2つ以上つけられるようになった
%% initiallize workspace
clear 
close all
%% Controlled node
c_n = {[13,30],[5,17],[1,3,10,24],[6,22,26,27],[14,15,16,29],[7,9,20,21,28]};
int_sys = cell(1,numel(c_n));
for i = 1 : numel(c_n)
    int_sys(i) = {strcat('controlled',num2str(i))};
end
int_sys = string(int_sys);
%% Mesurement Node
m_n = [1:30];%好きなだけ

%% Noise Port Node
n_n = [2];%1nodeのみ

%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
% 20 node 
% node parameter equal
% Cenect ratio 0.1
n = network_swing_simple(30, [1, 2], [2, 10]*1e-2, 1, [1,3], 0.1, 10);
n.Adj_ref = n.Adj_ref*0;

%% LQR Parameter
cell_Q = {diag([1, 1000, 1, 1000]),diag([1, 1000, 1, 1000]),diag([1, 1000, 1, 1000, 1, 1000, 1, 1000]),diag([1, 1000, 1, 1000, 1, 1000, 1, 1000]),diag([1, 1000, 1, 1000, 1, 1000, 1, 1000]),diag([1, 1000, 1, 1000, 1, 1000, 1, 1000, 1, 1000])};
cell_R = {diag([1e-3,1e-3]),diag([1e-3,1e-3]),diag([1e-3,1e-3,1e-3,1e-3]),diag([1e-3,1e-3,1e-3,1e-3]),diag([1e-3,1e-3,1e-3,1e-3]),diag([1e-3,1e-3,1e-3,1e-3,1e-3])};
%
%% add I/O port (Original System)
sys_org = n.get_sys();
for i =  1 : numel(c_n)
    idx = cell2mat(c_n(i));
    for j = 1 : numel(idx)
        sys_org = n.add_io(sys_org,idx(j), strcat(int_sys(i),'_node',num2str(idx(j))));
    end
end
%noise_port_idx = setdiff(cell2mat(n_n),cell2mat(c_n),'sorted');
for i = 1 :numel(n_n)
    sys_org = n.add_io(sys_org,n_n(i),strcat('node',num2str(n_n(i))));
end
%% Add Controller (Retro fit)
% Retro fit
for i = 1: numel(c_n)
    if iscell(c_n(i))
        idx = cell2mat(c_n(i));
    end
    Q = cell2mat(cell_Q(i));
    R = cell2mat(cell_R(i));
    n.add_controller(idx, Q, R);
end
sys_c1 = n.get_sys_controlled(sys_org);
%% Add Controller (Extend Retro fit)
% Obj.controllers initialize
n.controllers = {};
%全システム
sys_all = n.get_sys();
for i = 1 : numel(c_n)
    if iscell(c_n(i))
        idx = cell2mat(c_n(i));
    end
    %環境モデルの抽出
    [sys_local, sys_env] = n.get_sys_local(idx);
    model1 = balred(sys_env, 6);
    
    Q = cell2mat(cell_Q(i));
    R = cell2mat(cell_R(i));
    n.add_controller(idx,model1, Q, R);
    
end

sys_c2 = n.get_sys_controlled(sys_org);

%% Response Noise
%for i = 1 : 20
 %   n_n =i;
    close all
    Ts = 0.01;
    %Noise_Pattern Impulse or Colored noise
    N_P = 0;% CN:1% Im:0
    color = 'white';
    noise_power = 0.1;
    %Phase:1 or Fruquency:2
    check = 2;


    s_t = 2000;%simulation_time
    %noise = randn(2,s_t/Ts);
    time = 0:Ts:s_t-Ts;

    if N_P == 0
        % Insert Noise for Original System
        [y0, t0] = impulse(sys_org({'x'}, {strcat('d_node',num2str(n_n))}), time);
        % Insert Noise for retro fit controller of Node1
        [y1, t1] = impulse(sys_c1({'x', 'xhat_controlled1', 'u_controlled1'}, {strcat('d_node',num2str(n_n))}), time);
        % Insert Noise for Extend retro fit controller of Node1
        [y2, t2] = impulse(sys_c2({'x', 'xhat_controlled1', 'u_controlled1'}, {strcat('d_node',num2str(n_n))}), time);
    else
        % Generate noise
        cn = dsp.ColoredNoise('Color',color,'SamplesPerFrame',s_t/Ts,'NumChannels',1);
        noise = cn()*noise_power;
        noise = [zeros(s_t/Ts,1),noise];
        % Insert Noise for Original System
        [y0, t0] = lsim(sys_org({'x'}, {strcat('d_node',num2str(n_n))}), noise, time);
        % Insert Noise for retro fit controller of Node1
        [y1, t1] = lsim(sys_c1({'x', 'xhat_controlled1', 'u_controlled1'}, {strcat('d_node',num2str(n_n))}), noise, time);
        % Insert Noise for Extend retro fit controller of Node1
        [y2, t2] = lsim(sys_c2({'x', 'xhat_controlled1', 'u_controlled1'}, {strcat('d_node',num2str(n_n))}), noise, time);
    end

    % 描画範囲
    F_P = [600,300,600,600];
    limit_time = 300;
    l_t = limit_time/Ts;
    % c_n pointをすべてバラス
    c_n = cell2mat(c_n);    
    for i = 1 : numel(m_n)
        if (sum(double(n_n == m_n(i))) == 1) && (sum(double(c_n == m_n(i))) == 1)
            figure('position',F_P,'Name',sprintf('Response node %d:Noise Port & Controlled\n',m_n(i)))
        elseif sum(double(n_n == m_n(i))) == 1
            figure('position',F_P,'Name',sprintf('Response node %d:Noise Port\n',m_n(i)))   
        elseif sum(double(c_n == m_n(i))) == 1
            figure('position',F_P,'Name',sprintf('Response node %d:Controlled\n',m_n(i)))
        else
            figure('position',F_P,'Name',sprintf('Response node %d\n',m_n(i)))
        end
        Idx = [2*m_n(i)-1,2*m_n(i)];
        subplot(3,1,1)
        grid on
        plot(t0(2:l_t),y0(2:l_t,Idx,check),'LineWidth',1.5)
        legend('\theta','\omega','location','best')
        subplot(3,1,2)
        grid on
        plot(t1(2:l_t),y1(2:l_t,Idx,check),'LineWidth',1.5)
        subplot(3,1,3)
        grid on
        plot(t2(2:l_t),y2(2:l_t,Idx,check),'LineWidth',1.5)
    end
    
    %saveas(gcf,strcat('n_n_',num2str(n_n),'-c_n_',num2str(c_n),'.png'))
    %saveas(gcf,strcat('n_n_',num2str(n_n),'-c_n_',num2str(c_n),'.emf'))
    
%end

n.plot()
%saveas(gcf,'GraphStructure.png')