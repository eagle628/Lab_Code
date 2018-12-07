clear
close all

s = tf('s');

% G = (2*s+3)/(s+4);
% G = rss(2);
kkk = 6;
% G = 10*(s+5)*(s+50)/(s+10)/(s+100);
while true
    G = tf(10*[rand(1,kkk)],[1,rand(1,kkk-1)]);
    if isstable(G)
        break;
    end
end

% fbode(G);

N = 10000;
Ts = 0.01;

G_d = c2d(G, Ts, 'foh');
hhh = 5;
while true
    H_d = tf([1,randn(1,hhh)],[1,randn(1,hhh)],Ts);
    if isstable(H_d)
        break;
    end
end

% % aaa = order(G);
% % bbb = aaa+1;
% % 
% % ccc = order(H);
% % ddd = order(H);
% % 
% % init_sys = rss(aaa);
% % % lambda = 1e-6;
% % % init_sys = 1/(s+lambda)^aaa;
% % [init_num, init_den] = tfdata(init_sys,'v');
% % 
% % init_rho = [init_den(2:end),init_num];
% % rho = init_rho;

figure('position',[-900,0,600,600])
% fbode(G,[],':','linewidth',3.0);
bode(G)
hold on

condition = [kkk-1,kkk,hhh,hhh];

rng('shuffle')
seeds = randi(2^30,100,2);
idx = 1;
parfor_progress(100);
for itr = 1:100
    % input setting
    rng(seeds(itr,1))
    u = randn(N,1)*1;
    % noise setting
    rng(seeds(itr,2))
    d = randn(N,1)*1;
    y =  lsim(G_d, u) + lsim(H_d, d);
    data = iddata(y, u, Ts);
    
%     try
    sys = RIVC(data,condition);
%     fbode(sys.G,[]);
    bode(sys.G,'r');
%     catch ME
%         fprintf('%d th ERROR \n',idx);
%         idx = idx +1;
%     end
    drawnow
    parfor_progress();
end
parfor_progress(0);

% %  関数化前
% % % % % % % for itr1 = 1 : 5
% % % % % % %     den = [1, rho(1:aaa)];
% % % % % % %     num = rho(aaa+1:aaa+bbb);
% % % % % % % 
% % % % % % %     x = lsim(c2d(tf(num, den), Ts, 'foh'), u);
% % % % % % % 
% % % % % % %     y_fc = zeros(N, aaa+1);
% % % % % % %     x_fc = zeros(N, aaa+1);
% % % % % % %     for itr2 = 1 : aaa+1
% % % % % % %         num = zeros(size(num));
% % % % % % %         num(itr2) = 1;
% % % % % % %         y_fc(:, itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), y);
% % % % % % %         x_fc(:, itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), x);
% % % % % % %     end
% % % % % % %     u_fc = zeros(N, bbb);
% % % % % % %     for itr2 = 1 : bbb
% % % % % % %         num = zeros(size(num));
% % % % % % %         num(itr2) = 1;
% % % % % % %         u_fc(:, itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), u);
% % % % % % %     end
% % % % % % % 
% % % % % % %     xi = y-x;
% % % % % % %     noise_sys = armax(iddata(xi,[],Ts),[ccc,ddd]);
% % % % % % %     noise_sys = tf(noise_sys.A,noise_sys.C,Ts);
% % % % % % % %     noise_sys = tf(1,1,Ts);
% % % % % % % 
% % % % % % %     y_f = zeros(N, aaa+1);
% % % % % % %     x_f = zeros(N, aaa+1);
% % % % % % %     for itr2 = 1 : aaa+1
% % % % % % %         y_f(:, itr2) = lsim(noise_sys, y_fc(:, itr2));
% % % % % % %         x_f(:, itr2) = lsim(noise_sys, x_fc(:, itr2));
% % % % % % %     end
% % % % % % %     u_f = zeros(N, bbb);
% % % % % % %     for itr2 = 1 : bbb
% % % % % % %         u_f(:, itr2) = lsim(noise_sys, u_fc(:, itr2));
% % % % % % %     end
% % % % % % % 
% % % % % % %     phi = [-y_f(:,2:end), u_f];
% % % % % % % %     phi = [-x_f(:,2:end), u_f]; % noise free version
% % % % % % % 
% % % % % % %     pre_rho = (phi'*phi)\(phi'*y_f(:,1));
% % % % % % %     pre_rho = pre_rho';
% % % % % % % 
% % % % % % %     if norm(pre_rho-rho) < 1e-6 
% % % % % % %         break;
% % % % % % %     else
% % % % % % %         rho = pre_rho;
% % % % % % %     end
% % % % % % % 
% % % % % % % 
% % % % % % %     sys = tf(rho(aaa+1:aaa+bbb), [1,rho(1:aaa)]);
% % % % % % %     fbode(sys);
% % % % % % %     drawnow
% % % % % % % end
%% 書き直す前
% % % % % for itr1 = 1:10
% % % % %     den = [1,rho(1:aaa)];
% % % % %     num = zeros(size(den));
% % % % %     num(end-bbb+1:end) = rho(aaa+1:aaa+bbb);
% % % % %     % generate Insrumnetal varialbes
% % % % %     x = lsim(c2d(tf(num,den),Ts,'foh'), u);
% % % % %     % prefiltering u,y,& x
% % % % %     y_fc = zeros(N, aaa+1);
% % % % %     x_fc = zeros(N, aaa+1);
% % % % %     for itr2 = 1 : aaa+1
% % % % %         num = zeros(size(num));
% % % % %         num(end-itr2+1) = 1;
% % % % %         y_fc(:, end-itr2+1) = lsim(c2d(tf(num, den), Ts, 'foh'), y);
% % % % %         x_fc(:, end-itr2+1) = lsim(c2d(tf(num, den), Ts, 'foh'), x);
% % % % %     end
% % % % %     u_fc = zeros(N, bbb);
% % % % %     for itr2 = 1 : bbb
% % % % %         num = zeros(size(num));
% % % % %         num(end-itr2+1) = 1;
% % % % %         u_fc(:, end-itr2+1) = lsim(c2d(tf(num, den), Ts, 'foh'), u);
% % % % %     end
% % % % % 
% % % % %     % noise vector
% % % % %     xi = y-x;
% % % % % %     noise_sys = armax(iddata(xi,[],Ts),[ccc,ddd]);
% % % % % %     noise_sys = tf(noise_sys.A,noise_sys.C,Ts);
% % % % %     noise_sys = tf(1,1,Ts);
% % % % % 
% % % % %     y_f = zeros(N, aaa+1);
% % % % %     x_f = zeros(N, aaa+1);
% % % % %     for itr2 = 1 : aaa+1
% % % % %         num = zeros(size(num));
% % % % %         num(end-itr2+1) = 1;
% % % % %         y_f(:, end-itr2+1) = lsim(noise_sys, y_fc(:, end-itr2+1));
% % % % %         x_f(:, end-itr2+1) = lsim(noise_sys, x_fc(:, end-itr2+1));
% % % % %     end
% % % % %     u_f = zeros(N, bbb);
% % % % %     for itr2 = 1 : bbb
% % % % %         num = zeros(size(num));
% % % % %         num(end-itr2+1) = 1;
% % % % %         u_f(:, end-itr2+1) = lsim(noise_sys, u_fc(:, end-itr2+1));
% % % % %     end
% % % % % 
% % % % %     phi = [-y_f(:,2:end), u_f];
% % % % % 
% % % % %     pre_rho = (phi'*phi)\(phi'*y_f(:,1));
% % % % %     pre_rho = pre_rho';
% % % % % 
% % % % %     if norm(pre_rho-rho) < 1e-6 && itr1 > 5
% % % % %         break;
% % % % %     else
% % % % %         rho = pre_rho;
% % % % %     end
% % % % %     
% % % % %     sys = tf(rho(aaa+1:aaa+bbb), [1,rho(1:aaa)]);
% % % % %     fbode(sys);
% % % % %     drawnow
% % % % % end


% % function sys = SRIVC()
% % 
% % end
