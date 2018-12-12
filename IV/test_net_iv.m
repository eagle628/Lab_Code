clear
close all
%% genereate Network
seed = 10;
Node_number = 10;
n_ori = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% control & noise Node
c_n = 1;
n_n = [2:Node_number];
%% signal power
% id_in_p = 1;
% noise_p = 10;
%% Node Name
ob_v_p  = {strcat('v_node',num2str(c_n))};
ob_w_p  = {strcat('w_node',num2str(c_n))};
ob_y_p  = {strcat('y_node',num2str(c_n))};
ID_in_p = {strcat('d_node',num2str(c_n))};
ob_xhat_p = {'xhat_controlled1'}; % The last number should be Controllers Group number.
Noise = cell(1,numel(n_n));
for i = 1 : numel(n_n)
Noise(i) = {strcat('d_node',num2str(n_n(i)))};
end
%% add I/O port for identificaiton
sys_ori = n_ori.get_sys();
add_io_node = unique([c_n,n_n]);
for nnn = add_io_node
    sys_ori = n_ori.add_io(sys_ori, nnn, strcat('node',num2str(nnn)));
end
[sys_local, sys_env] = n_ori.get_sys_local(c_n);
sys_local_vw = sys_local({'w'},{'v'});
%% env character
[mag_ori,~,wout_ori] = bode(sys_env);
wrange = {wout_ori(1),wout_ori(end)};

%% vw generate
rng('shuffle')
id_in_p = 10;
noise_p = 1;

N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;


model_dim = 6;

max_itr = 100;
parfor_progress(max_itr);

mag_result_set = cell(max_itr, 1);
wout_result_set = cell(max_itr,1);

sys_v = sys_ori(ob_v_p, cat(2,ID_in_p,Noise));
sys_w = sys_ori(ob_w_p, cat(2,ID_in_p,Noise));
parfor itr = 1 : max_itr
    d = randn(N,2+numel(n_n)*2);
    d(:,1:2) = d(:,1:2)*id_in_p;
    d(:,3:end) = d(:,3:end)*noise_p;

    v = lsim(sys_v, d, t);
    w = lsim(sys_w, d, t);

    data = struct();
    data.Ts = Ts;
    data.u = w;
    data.y = v;

    try
%     G2 = CLIVC2(data,-sys_local_vw,[5,6],1);
    G2 = CLRIVC(data,-sys_local_vw,[model_dim-1,model_dim,model_dim-1,model_dim-1],100);
    [mag_result_set{itr},~,wout_result_set{itr}] = bode(G2, wrange);
    end
    parfor_progress();
end
parfor_progress(0);

figure
for itr = 1 : max_itr
    try
    semilogx(wout_result_set{itr},mag2db(squeeze(mag_result_set{itr})),'r');
    end
    hold on, grid on ,box  on
    drawnow
end
semilogx(wout_ori,mag2db(squeeze(mag_ori)),'b:','LineWidth',3.0);
ax = gca;
ax.XScale ='log';
