clear
close all
%% genereate Network
seed = 9;
Node_number = 4;
n_ori = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
% n_ori = network_swing_simple(Node_number, 1, [2,10]*1e-2, 1, 1, 0.1, seed);
n_ori.Adj_ref = n_ori.Adj_ref*0;
n_ori.plot()
%% control & noise Node
c_n = 1;
n_n = [2:Node_number];
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
Q = diag([1,1000]);
R = 1e-3;
% n_ori.add_controller(c_n, Q, R);
n_ori.add_controller(n_n);
sys_con = n_ori.get_sys_controlled(sys_ori);

% sys_con = sys_ori;

Loop = loopsens(sys_env,-sys_local_vw)
isstable(sys_con({'w_node1'},sys_ori.InputGroup.d_node1))
isstable(sys_con({'v_node1'},sys_ori.InputGroup.d_node1))
isstable(sys_con)



rng('shuffle')

id_p = 1;
noise_p = 0;

N = 1000;
Ts = 0.01;
t = (0:N-1)'*Ts;


model_dim = 6;

max_itr = 100;

parfor_progress(max_itr);

mag_result_set = cell(max_itr, 1);
wout_result_set = cell(max_itr,1);
data_set = cell(max_itr,1);
IDsys_set = cell(max_itr,1);

sys_v_id = sys_con(ob_v_p, ID_in_p);
sys_v_np = sys_con(ob_v_p, Noise);
sys_w_id = sys_con(ob_w_p, ID_in_p);
sys_w_np = sys_con(ob_w_p, Noise);
error = 0;
Error_message = {};
parfor itr = 1 : max_itr
    d_id = zeros(N, 2);
%     d_id = randn(N, 2);
    d_id(:,1) = randn(N, 1)*id_p;
    d_np = zeros(N,numel(n_n)*2);
    d_np(:,1:2:end) = randn(N,numel(n_n))*noise_p;

    v_id = lsim(sys_v_id, d_id, t);
    w_id = lsim(sys_w_id, d_id, t);

    v_np = lsim(sys_v_np, d_np, t);
    w_np = lsim(sys_w_np, d_np, t);
    
    v = v_id + v_np;
    w = w_id + w_np;
    
    data = struct();
    data.Ts = Ts;
    data.u = w;
    data.y = v;
    data.r = w - lsim(c2d(sys_local_vw, Ts, 'foh'), v);
    data.io = {v_id,v_np,w_id,w_np};
    data_set{itr} = data;

    try
%         G2 = CLIVC2(data,-sys_local_vw,[5,6],1);
        sys = CLRIVC(data,-sys_local_vw,[model_dim-1,model_dim,model_dim-1,model_dim-1],100);
%         sys = RIVC(data,[model_dim-1,model_dim,model_dim-1,model_dim-1],100);
        [mag_result_set{itr},~,wout_result_set{itr}] = bode(sys.G, wrange);
        IDsys_set{itr} = sys;
    catch ME
        error = error+1;
        Error_message = [Error_message,{ME}];
    end
    parfor_progress();
end
parfor_progress(0);

fprintf('Error number is %d.\n',error)

%%
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

%% 
rsme_set = zeros(1,max_itr);

for itr = 1 : max_itr
    try
        G_test = c2d(IDsys_set{itr}.G, Ts, 'foh');
        cloop_d = loopsens(G_test,-c2d(sys_local_vw,Ts));
        y_test = lsim(G_test*cloop_d.Si, data_set{itr}.r);
        rsme_set(itr) = norm(data_set{itr}.io{1} - y_test);
    catch
        rsme_set(itr) = nan;
    end
end
rsme_set = rsme_set(~isnan(rsme_set));

figure;
boxplot(rsme_set);
%% 
% BFR_set = zeros(1,max_itr);
% for itr = 1 : max_itr
%     try
%         G_test = c2d(IDsys_set{itr}.G, Ts, 'foh');
%         cloop_d = loopsens(G_test,-c2d(sys_local_vw,Ts));
%         y_test = lsim(G_test*cloop_d.Si, data_set{itr}.r);
%         BFR_set(itr) = 1-norm(data_set{itr}.io{1} - y_test)/norm(data_set{itr}.io{1} - mean(data_set{itr}.io{1}));
%     catch
%         BFR_set(itr) = nan;
%     end
% end
% 
% figure;
% boxplot(BFR_set);
%%
t = (0:1:N-1)'*Ts;
figure('name','I/O data');
plot(t,data_set{1}.u,'b',t,data_set{1}.y,'r')
legend('Input','Output')
figure('Name','output S/N');
plot(t,data_set{1}.io{1},'b',t,data_set{1}.io{2},'g')
legend('Output : S','Output : N')
