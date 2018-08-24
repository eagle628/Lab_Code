%{
このコードはレトロフィット制御器の相互信号を雑音としてサブシステムに認識させれば
環境にノイズを入れた場合でも問題何のではないかというコンセプトで行う．
相互信号が計測可能であるとの前提の下で行う．

無理気な気配がぷんぷん
無理でした．　2018/05/21

%}
%% initiallize workspace
clear 
close all
%% Controlled node
c_n = {[1]};
%c_n = {[13,30],[5,17],[1,3,10,24],[6,22,26,27],[14,15,16,29]};
int_sys = cell(1,numel(c_n));
for i = 1 : numel(c_n)
    int_sys(i) = {strcat('controlled',num2str(i))};
end
int_sys = string(int_sys);
%% Mesurement Node
m_n = [1];%好きなだけ

%% Noise Port Node
n_n = [3];%1nodeのみ

%% Network structure
% Randam Number Control (Seed)
rng(10);
% Generate Network (x,u,v)
% 20 node 
% node parameter equal
% Cenect ratio 0.1
n = network_swing_simple(30, [1, 2], [2, 10]*1e-2, 1, [1,5], 0.1, 10);
n.Adj_ref = n.Adj_ref*0;

%% LQR Parameter
cell_Q = {diag([1, 1000])};
cell_R = {diag([1e-3])};
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
Retro_controller = n.controllers;
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
    close all
    Ts = 0.01;
    %Noise_Pattern Impulse or Colored noise
    N_P = 0;% CN:1% Im:0
    color = 'white';
    noise_power = 0.1;
    %Phase:1 or Fruquency:2
    check = 1;


    s_t = 300;%simulation_time
    %noise = randn(2,s_t/Ts);
    time = 0:Ts:s_t-Ts;

    input = zeros(2,length(time));
    input(2,1) = 100;
    [y0, t0] = lsim(sys_org({'x'}, {strcat('d_node',num2str(n_n))}), input, time, 'zoh');
    [y0_2, t0_2] = lsim(sys_c1({'x','v_controlled1'}, {strcat('d_node',num2str(n_n))}),input ,time, 'zoh');
    
    
    sys_converse = n.sys_inter_to_dis(sys_c1,cell2mat(c_n(1)));
    [y0_3, t0_3] = lsim(sys_c1({'x','v_controlled1'}, {strcat('d_node',num2str(n_n))}), input, time,'zoh');
    [y0_4, t0_4] = lsim(sys_converse({'x','v_controlled1'}, {strcat('d_node',num2str(n_n))}), input, time,'zoh');
    
    %{
    y0_5 = zeros(length(time),63);
    t0_5 = zeros(length(time),1);
    for i = 1:(length(time)-1)
        if i == 1
            y0_5([i,i+1],:) = lsim(sys_c1({'x', 'xhat_controlled1','v_controlled1'}, {strcat('d_node',num2str(n_n)),'d_controlled1'}),[100,0;0,0;0,0;0,0;],[0,Ts],'zoh');
        else
            %y0_5([i,i+1],:) = lsim(sys_c1({'x', 'xhat_controlled1','v_controlled1'}, {strcat('d_node',num2str(n_n)),'d_controlled1'}),[0,0;0,0;0,y0_5([i],63);0,0;],[0,Ts],y0_5([i],1:62),'zoh');
            y0_5([i,i+1],:) = lsim(sys_c1({'x', 'xhat_controlled1','v_controlled1'}, {strcat('d_node',num2str(n_n)),'d_controlled1'}),[0,0;0,0;0,0;0,0;],[0,Ts],y0_5([i],1:62),'zoh');
        end
    end
    %}
    %%
    %{
    i = 1;
    figure
    subplot(3,1,1)
    plot(t0,y0_3(:,[2*i-1,2*i],1))
    xlim([0 300])
    subplot(3,1,2)
    plot(t0,y0_4(:,[2*i-1,2*i]))
    xlim([0 300])
    subplot(3,1,3)
    plot(t0,y0_5(:,[2*i-1,2*i]))
    xlim([0 300])
    %}
    %
    for i = 1
        figure
        subplot(3,1,1)
        plot(t0,y0(:,[2*i-1,2*i],1))
        xlim([0 300])
        subplot(3,1,2)
        plot(t0,y0_3(:,[2*i-1,2*i]))
        xlim([0 300])
        subplot(3,1,3)
        plot(t0,y0_4(:,[2*i-1,2*i]))
        xlim([0 300])
    end
    %}
    %{
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
        noise = [noise,noise];
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
    
    %}