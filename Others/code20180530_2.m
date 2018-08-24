% last version
%% initiallize workspace
clear 
close all
%% Controlled node
%c_n = {[2,25,26,27,28,14,30]};
%c_n = {[3,26,30,13]};
%c_n = {[15,16,22,24,29]};
%c_n = {[3,26,30,13],[15,16,22,24,29],[29,30]};
%c_n = {[2,3,9],[1,8],[4,7],[5,6,10]};
%c_n = {[1:10],[11,20],[21:30]};
c_n = {1};
int_sys = cell(1,numel(c_n));
for i = 1 : numel(c_n)
    int_sys(i) = {strcat('controlled',num2str(i))};
end
int_sys = string(int_sys);
%% Mesurement Node
m_n = [1];%好きなだけ

%% Noise Port Node
n_n = [1];%1nodeのみ

%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
% 20 node 
% node parameter equal
% Cenect ratio 0.1
n = network_swing_simple(10, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, 10);
n.Adj_ref = n.Adj_ref*0;

%% LQR Parameter
cell_Q = cell(1,numel(c_n));
cell_R = cell(1,numel(c_n));
for i =  1 : numel(c_n)
    s = length(cell2mat(c_n(i)));
    cell_Q(i) = {kron(eye(s),diag([1,1000]))};
    cell_R(i) = {kron(eye(s),diag([1e-3]))};
%cell_Q = {diag([1, 1000, 1, 1000, 1, 1000, 1, 1000, 1, 1000]),diag([1, 1000, 1, 1000])};
%cell_R = {diag([1e-3,1e-3,1e-3,1e-3,1e-3]),diag([1e-3,1e-3])};
end
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
% Obj.controllers initialize
n.controllers = {};
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

%% Set Noise
%for i = 1 : 20
 %   n_n =i;
    close all
    Ts = 0.01;
    %Noise_Pattern Impulse or Colored noise
    N_P = 1;% CN:1% Im:0 % Cos:2
    % Color Noise
    color = 'white';
    noise_power = 0.1;
    % Cos Fruquency 
    noise_omega = 10^2;
    
    %Phase:1 or Fruquency:2
    check = 2;
    
%% Calculate Respomse

    s_t = 2000;%simulation_time
    %noise = randn(2,s_t/Ts);
    time = 0:Ts:s_t-Ts;

    %noise port
    noise_port = cell(1,numel(n_n));
    for i = 1: numel(n_n)
        noise_port(i) = {strcat('d_node',num2str(n_n(i)))};
    end
  
    if N_P == 0
        % Insert Noise for Original System
        [y0, t0] = impulse(sys_org({'x'}, noise_port), time);
        % Insert Noise for retro fit controller of Node1
        [y1, t1] = impulse(sys_c1({'x', 'xhat_controlled1', 'u_controlled1'}, noise_port), time);
        % Insert Noise for Extend retro fit controller of Node1
        [y2, t2] = impulse(sys_c2({'x', 'xhat_controlled1', 'u_controlled1'}, noise_port), time);
    else
        % Generate noise
        if N_P == 1
            cn = dsp.ColoredNoise('Color',color,'SamplesPerFrame',s_t/Ts,'NumChannels',1);
            noise = cn()*noise_power;
        else
            noise = cos(2*pi/noise_omega.*time);
            noise = noise';
        end
        if check == 2
            noise = [zeros(s_t/Ts,1),noise];
        else
            noise = [noise,zeros(s_t/Ts,1)];
        end
        % Insert Noise for Original System
        [y0, t0] = lsim(sys_org({'x'}, {strcat('d_node',num2str(n_n))}), noise, time);
        % Insert Noise for retro fit controller of Node1
        [y1, t1] = lsim(sys_c1({'x', 'xhat_controlled1', 'u_controlled1'}, {strcat('d_node',num2str(n_n))}), noise, time);
        % Insert Noise for Extend retro fit controller of Node1
        [y2, t2] = lsim(sys_c2({'x', 'xhat_controlled1', 'u_controlled1'}, {strcat('d_node',num2str(n_n))}), noise, time);
    end

    %% plot area
    % 描画範囲
    F_P = [600,300,600,600];
    limit_time = 300;
    l_t = limit_time/Ts;
    % c_n pointをすべてバラス
    c_n = cell2mat(c_n);    
    for i = 1 : numel(m_n)
        if (sum(double(n_n == m_n(i))) >= 1) && (sum(double(c_n == m_n(i))) >= 1)
            figure('position',F_P,'Name',sprintf('Response node %d:Noise Port & Controlled\n',m_n(i)))
        elseif sum(double(n_n == m_n(i))) >= 1
            figure('position',F_P,'Name',sprintf('Response node %d:Noise Port\n',m_n(i)))   
        elseif sum(double(c_n == m_n(i))) >= 1
            figure('position',F_P,'Name',sprintf('Response node %d:Controlled\n',m_n(i)))
        else
            figure('position',F_P,'Name',sprintf('Response node %d\n',m_n(i)))
        end
        Idx = [2*m_n(i)-1,2*m_n(i)];
        if size(y0,3) == 1
            check = 1;
        end
        subplot(3,1,1)
        grid on;hold on;
        y = y0(1:l_t,Idx,check);
        %y = sum(y0(1:l_t,Idx,[2,4]),3);
        plot(t0(1:l_t),y,'LineWidth',1.5)
        legend('\theta','\omega','location','best')
        subplot(3,1,2)
        grid on;hold on;
        y = y1(1:l_t,Idx,check);
        %y = sum(y1(1:l_t,Idx,[2,4]),3);
        plot(t1(1:l_t),y,'LineWidth',1.5)
        subplot(3,1,3)
        grid on;hold on;
        y = y2(1:l_t,Idx,check);
        %y = sum(y2(1:l_t,Idx,[2,4]),3);
        plot(t2(1:l_t),y,'LineWidth',1.5)
    end
    
    %saveas(gcf,strcat('n_n_',num2str(n_n),'-c_n_',num2str(c_n),'.png'))
    %saveas(gcf,strcat('n_n_',num2str(n_n),'-c_n_',num2str(c_n),'.emf'))
    
%end

%n.plot()
%saveas(gcf,'GraphStructure.png')