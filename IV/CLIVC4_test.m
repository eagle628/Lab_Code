%% NOTE
% Identification of Continuous-time Models from Sampled Data.pdf 
% Chapter 5
% section 5.5.3
% Four Step CLIVC4 algorithm
%% main
clear
close all
%% Prepare Plan and Controller and Noise Filter
% 
Flag1 = true;

if Flag1 
    rng(20)
    GGG = 7;
    while true
%         G = tf([5*rand(1,GGG-1),1e-15],[1,10*rand(1,GGG-1)]);
        G = tf(5*[randn(1,GGG-1),0],[1,10*rand(1,GGG-1)]);
%         G = tf([-5*rand(1),5*randn(1,GGG-2),1e-10*randn(1)],[1,10*rand(1,GGG-1)]);
%         G = tf([-5*rand(1),5*randn(1,GGG-2),0],[1,10*rand(1,GGG-1)]);
    %     G = rss(GGG-1);
    %     G.D = randn(1);
%         if isstable(G) && (sum(real(zero(G)>0)) == 0)
%         if (sum(real(zero(G)>0)) == 1)
%         G = rss(GGG-1);
%         [~,~,~,d] = ssdata(G);
%         if d ~= 0
            if isstable(G)
                break;
            end
%         end
    end
    disp('G_Zero')
    zero(G)
    disp('G_Pole')
    pole(G)

    % % load ex_network_sys5.mat sys_env sys_local_vw
    % % G = balreal(sys_env;

    rng(28)
    CCC = 3;
    while true
%         C = tf([rand(1,1)],[1,rand(1,CCC-2),0]);
        C = rss(CCC);
        feed = loopsens(G, C);
%         if ~feed.Stable && (sum(feed.Poles>0) == 1) && (max(real(feed.Poles))<1e-8)
        if feed.Stable
            disp('Loop_Pole')
            feed.Poles
            break;
        end
    end
    % C = -sys_local_vw;

    N = 10000;
    Ts = 0.01;

    G_d = c2d(G, Ts, 'foh');
    C_d = c2d(C, Ts, 'foh');

    rng(1024)
    HHH = GGG-1;
    while true
        H_d = tf([1,randn(1,HHH)],[1,randn(1,HHH)],Ts);
        if isstable(H_d)
            break;
        end
    end
%     H = ss(G);
    % H.B = zeros(GGG-1, 1);
    % H.B(4) = 1;

%     H_d = c2d(H, Ts, 'foh');

else
    load 
end

na = order(G);
nb = order(G)+1;
nc = na;
nd = na;
%% Generate Response : Availabel r1 r2 y u
disp('Generate I/O data')

max_itr = 100;
data_set = cell(1,max_itr);
rng('shuffle')

parfor_progress(max_itr);
parfor itr = 1 : max_itr
    r1 = randn(N, 1)*0;
    r2 = randn(N, 1)*1;
    d  = randn(N, 1)*0;
    r  = r1 + lsim(C_d, r2);
    cloop_d = loopsens(G_d,C_d);
    yyy = lsim(G_d*cloop_d.Si, r);
    ddd = lsim(cloop_d.So*H_d, d);
    y = yyy + ddd;
    u = r - lsim(C_d, y);
    data = struct();
    data.u = u;
    data.r = r;
    data.y = y;
    data.yyy = yyy;
    data.ddd = ddd;
    data.Ts = Ts;
    data_set{itr} = data;
    parfor_progress();
end
parfor_progress(0);
%%
disp('IV Method')

[mag_ori,~,wout_ori] = bode(G);
% wrange = {wout_ori(1),wout_ori(end)};
wrange = wout_ori;
mag_result_set = cell(1,max_itr);
wout_result_set = cell(1,max_itr);
IDsys_set = cell(1,max_itr);
parfor_progress(max_itr);
error = 0;
parfor itr = 1 : max_itr
    try
%         IDsys_set{itr} = CLIVC2(data_set{itr},C_d,[na,nb],1);
%         [mag,~,wout] = bode(IDsys_set{itr},wrange);
        IDsys_set{itr} = CLRIVC(data_set{itr},C_d,[na,nb,nc,nd],1,10);
%         [num,den] = tfdata(G);
%         IDsys_set{itr} = CLRIVC_opt(data_set{itr},C_d,[na,nb,nc,nd],inv(tf(den,1)*d2c(H_d)));
        [mag,~,wout] = bode(IDsys_set{itr}.G,wrange);
        mag_result_set{itr} = mag;
        wout_result_set{itr} = wout;
    catch ME
        error = error + 1;
    end
    parfor_progress();
end
parfor_progress(0);

fprintf('Number of Error is %d.\n',error)
%% Drawing
figure
for itr = 1 : max_itr
    try
        semilogx(wout_result_set{itr}, mag2db(squeeze(mag_result_set{itr})),'r-');
        hold on, box on, grid on
        drawnow
    end
end
semilogx(wout_ori, mag2db(squeeze(mag_ori)),'b:','LineWidth',3.0);
ax = gca;
ax.XScale ='log';

%% 
% Identification of dynamic networks operating in the presence of algebraic loops.pdf
% BFR

rsme_set = zeros(1,max_itr);

for itr = 1 : max_itr
    try
        G_test = c2d(IDsys_set{itr}.G, Ts, 'foh');
        cloop_d = loopsens(G_test,C_d);
        y_test = lsim(G_test*cloop_d.Si, data_set{itr}.r);
        rsme_set(itr) = norm(data_set{itr}.yyy - y_test);
    catch
        rsme_set(itr) = inf;
    end
end
rsme_set = rsme_set(~isinf(rsme_set));

figure;
boxplot(rsme_set);
%%
t = (0:1:N-1)'*Ts;
figure('name','I/O data');
plot(t,data_set{1}.u,'b',t,data_set{1}.y,'r')
legend('Input','Output')
figure('Name','output S/N');
plot(t,data_set{1}.yyy,'r',t,data_set{1}.ddd,'g')
legend('Output : S', 'Output : N ')
