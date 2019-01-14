clear
close all

% % % table = struct();
% % % table(1).seed = 1;
% % % table(1).c_n = [1,2];
% % % table(2).seed = 2;
% % % table(2).c_n = [1,3];
% % % table(3).seed = 4;
% % % table(3).c_n = [1,2,4];
% % % table(4).seed = 8;
% % % table(4).c_n = 1;
% % % 
% % % base_dir = 'C:\Users\NaoyaInoue\Desktop\figure_set\CLIV_Method\various';
% % % for itr1 = 1 : length(table)
% % %     for itr2 = 1 : length(table(itr1).c_n)
% % %         for stable = [0, 1]
% % %             name = strcat(...
% % %                         'data_seed',num2str(table(itr1).seed),...
% % %                         '_c_n',num2str(table(itr1).c_n(itr2)),...
% % %                         '_stable',num2str(stable)...
% % %                         );
% % %             location = strcat(base_dir, '\', name);
% % %             iter_func(table(itr1).seed, table(itr1).c_n(itr2), stable, location);
% % %         end
% % %     end
% % % end
% % % 
% % % function iter_func(seed, c_n, Flag1, location)
%% genereate Network
seed = 1;
Node_number = 4;
net1 = network_swing_simple(Node_number, [1,2], [2,10]*1e-2, 1, [1,5], 0.1, seed);
% net1 = network_swing_simple(Node_number, 1, [2,10]*1e-2, 1, [0,1], 0.1, seed);

%% control & noise Node
c_n = 2;
n_n = [];
% 
% tmp_idx = find(net1.Adj(c_n, :)~=0);
% tmp_adj = 2 + 5*rand(1, length(tmp_idx));
% net1.Adj(c_n ,tmp_idx) = tmp_adj;
% net1.Adj(tmp_idx, c_n) = tmp_adj;
% net1.Adj_ref = net1.Adj_ref*0;
net1.plot()
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
sys_all = net1.get_sys();
add_io_node = unique([c_n,n_n]);
for nnn = add_io_node
    sys_all = net1.add_io(sys_all, nnn, strcat('node',num2str(nnn)));
end

Flag1 = false;

if Flag1
    % sys_stable = sys_all([1:8],[2:4]);
    not_idx = setdiff(1:2*net1.N, [2*c_n,2*c_n-1]);
    sys_stable = sys_all({'x'},{'u'});
    K = lqr(sys_stable, eye(size(sys_stable.A,1)), eye(size(sys_stable.B,2)));
    % K = lqr(sys_stable, eye(8), eye(3));
    K = sys_stable.B*K;
    K = blkdiag(zeros(2), K(not_idx, not_idx)); % environment adjustment
    sys_stable = sys_all;
    sys_stable.A  = sys_stable.A - K;

    [sys_local, sys_env] = net1.get_sys_local(c_n);
    sys_env.A = sys_env.A - K(not_idx, not_idx); % environment adjustment
    sys_local_vw = sys_local({'w'},{'v'});
else
    sys_stable = sys_all;
    [sys_local, sys_env] = net1.get_sys_local(c_n);
    sys_local_vw = sys_local({'w'},{'v'});
end
%% env character
[mag_ori,~,wout_ori] = bode(sys_env);
wrange = {wout_ori(1),wout_ori(end)};

Loop = loopsens(sys_env,-sys_local_vw);
Loop.Stable
isstable(sys_stable(ob_w_p,ID_in_p))
isstable(sys_stable(ob_v_p,ID_in_p))
isstable(sys_stable)
%% vw generate
Q = diag([1,1000]);
R = 1e-3;
% net1.add_controller(c_n, Q, R);
% net1.add_controller(n_n);
% net1.add_controller(c_n);
sys_con = net1.get_sys_controlled(sys_stable);


rng('shuffle')

id_p = 1;
noise_p = 0;

N = 10000;
Ts = 0.01;
t = (0:N-1)'*Ts;


% model_dim = 2*(Node_number-1);
model_dim = 6;

max_itr = 100;

parfor_progress(max_itr);

data_set = cell(max_itr,1);

iv_mag_result_set = cell(max_itr, 1);
iv_wout_result_set = cell(max_itr,1);
iv_IDsys_set = cell(max_itr,1);

oe_mag_result_set = cell(max_itr, 1);
oe_wout_result_set = cell(max_itr,1);
oe_IDsys_set = cell(max_itr,1);

% sys_v_id = balreal(sys_con(ob_v_p, ID_in_p));
% sys_w_id = balreal(sys_con(ob_w_p, ID_in_p));
% sys_v_id = sys_con(ob_v_p, ID_in_p);
% sys_w_id = sys_con(ob_w_p, ID_in_p);
[sys_v_id, ~, Tv, Tiv] = balreal(c2d(sys_con(ob_v_p, ID_in_p), Ts, 'foh'));
[sys_w_id, ~, Tw, Tiw] = balreal(c2d(sys_con(ob_w_p, ID_in_p), Ts, 'foh'));



Flag2 = isempty(Noise);
if ~Flag2
%     sys_v_np = balreal(sys_con(ob_v_p, Noise));
%     sys_w_np = balreal(sys_con(ob_w_p, Noise));
    sys_v_np = sys_con(ob_v_p, Noise);
    sys_w_np = sys_con(ob_w_p, Noise);
else
    sys_v_np = ss([],[],[],zeros(1, 2));
    sys_w_np = ss([],[],[],zeros(1, 2));
end
error = 0;
Error_message = {};
parfor itr = 1 : max_itr
%     d_id = zeros(N, 2);
%     d_id(:, 2) = randn(N, 1)*id_p;
    d_id = randn(N, 2);

%     d_np = randn(N,numel(n_n)*2)*noise_p;
    d_np = zeros(N,numel(n_n)*2);
%     d_np(:,2:2:end) = randn(N,numel(n_n))*noise_p;

    v_id = lsim(sys_v_id, d_id, t);
    w_id = lsim(sys_w_id, d_id, t);

    if ~Flag2
        v_np = lsim(sys_v_np, d_np, t);
        w_np = lsim(sys_w_np, d_np, t);
    else
        v_np = zeros(size(v_id));
        w_np = zeros(size(w_id));
    end
    
    v = v_id + v_np;
    w = w_id + w_np;
    
    data = struct();
    data.Ts = Ts;
%     data.u = movmean(w, 5);
%     data.y = movmean(v, 5);
%     data.r = movmean(w - lsim(c2d(sys_local_vw, Ts, 'foh'), v), 5);
    data.u = w;
    data.y = v;
    data.r = w - lsim(c2d(sys_local_vw, Ts, 'foh'), v);
    data.io = {v_id,v_np,w_id,w_np};
    data_set{itr} = data;

    try
%         G2 = CLIVC2(data,-sys_local_vw,[5,6],1);
%         sys = CLRIVC(data,balreal(-sys_local_vw),[model_dim,model_dim+1,model_dim,model_dim],1,10);
        sys = CLRIVC(data,balreal(-sys_local_vw),[model_dim,model_dim+1,model_dim,model_dim],1,10);
%         sys = SPEM_CLRIVC(data,-sys_local_vw,[model_dim,model_dim+1,model_dim,model_dim],1);
%         sys = RIVC(data,[model_dim-1,model_dim,model_dim-1,model_dim-1],100);
        [iv_mag_result_set{itr},~,iv_wout_result_set{itr}] = bode(sys.G, wout_ori);
        iv_IDsys_set{itr} = sys;
        oe_IDsys_set{itr} = oe(iddata(v, w, Ts), [model_dim+1, model_dim, 0]);
        [oe_mag_result_set{itr},~,oe_wout_result_set{itr}] = bode(oe_IDsys_set{itr}, wout_ori);
    catch ME
        error = error+1;
        Error_message = [Error_message,{ME}];
    end
    parfor_progress();
end
parfor_progress(0);

fprintf('Error number is %d.\n',error)
% % % save(location)
% % % end

%%
figure
% % 
% % for itr = 1 : max_itr
% %     try
% %     semilogx(iv_wout_result_set{itr},mag2db(squeeze(iv_mag_result_set{itr})),'r');
% %     semilogx(oe_wout_result_set{itr},mag2db(squeeze(oe_mag_result_set{itr})),'g');
% %     end
% %     hold on, grid on ,box  on
% %     drawnow
% % end

iv_result = [];
oe_result = [];
for itr = 1 : max_itr
    try
    iv_result = [iv_result,squeeze(iv_mag_result_set{itr})];
    oe_result = [oe_result,squeeze(oe_mag_result_set{itr})];
    end
end
hold on, grid on ,box  on
semilogx(wout_ori, mag2db(mean(iv_result,2)),'r')
semilogx(wout_ori, mag2db(mean(oe_result,2)),'g')

semilogx(wout_ori,mag2db(squeeze(mag_ori)),'b:','LineWidth',3.0);
ax = gca;
ax.XScale ='log';

legend('CL-RIVC','OE','Original','location','best')

%% 
iv_rsme_set = zeros(1,max_itr);
oe_rsme_set = zeros(1,max_itr);

for itr = 1 : max_itr
    try
        G_test = c2d(iv_IDsys_set{itr}.G, Ts, 'foh');
        cloop_d = loopsens(G_test,-c2d(sys_local_vw,Ts));
        y_test = lsim(G_test*cloop_d.Si, data_set{itr}.r);
        iv_rsme_set(itr) = norm(data_set{itr}.io{1} - y_test);
        G_test = tf(oe_IDsys_set{itr});
        cloop_d = loopsens(G_test,-c2d(sys_local_vw,Ts));
        y_test = lsim(G_test*cloop_d.Si, data_set{itr}.r);
        oe_rsme_set(itr) = norm(data_set{itr}.io{1} - y_test);
    catch ME
        iv_rsme_set(itr) = inf;
        oe_rsme_set(itr) = inf;
    end
end
iv_rsme_set = iv_rsme_set(~isinf(iv_rsme_set));
oe_rsme_set = oe_rsme_set(~isinf(oe_rsme_set));

figure;
boxplot([iv_rsme_set',oe_rsme_set'],'Labels',{'iv','oe'});
%%
t = (0:1:N-1)'*Ts;
figure('name','I/O data');
plot(t,data_set{1}.u,'b',t,data_set{1}.y,'r')
legend('Input','Output')
figure('Name','output S/N');
plot(t,data_set{1}.io{1},'b',t,data_set{1}.io{2},'g')
legend('Output : S','Output : N')
