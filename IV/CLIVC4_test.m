%% NOTE
% Identification of Continuous-time Models from Sampled Data.pdf 
% Chapter 5
% section 5.5.3
% Four Step CLIVC4 algorithm
%% main
clear
close all
%% Prepare Plan and Controller and Noise Filter
rng(20)
GGG = 4;
while true
    G = tf([rand(1,GGG)],[1,rand(1,GGG-1)]);
    if isstable(G)
        break;
    end
end

na = GGG-1;
nb = GGG;

rng(6)
CCC = 2;
while true
    C = tf([rand(1,CCC)],[1,rand(1,CCC-1)]);
    feed = loopsens(G, C);
    if isstable(feed)
        break;
    end
end

N = 10000;
Ts = 0.01;

G_d = c2d(G, Ts, 'foh');
C_d = c2d(C, Ts, 'foh');

rng(1024)
HHH = 4;
while true
    H_d = tf([1,randn(1,HHH)],[1,randn(1,HHH)],Ts);
    if isstable(H_d)
        break;
    end
end

nc = HHH;
nd = HHH;

%% Generate Response : Availabel r1 r2 y u
max_itr = 20;
data_set = cell(1,max_itr);
rng('shuffle')

parfor_progress(max_itr);
parfor itr = 1 : max_itr
    r1 = randn(N, 1)*0;
    r2 = randn(N, 1)*1;
    d  = randn(N, 1)*0.1;
    r  = r1 + lsim(C_d, r2);
    cloop_d = loopsens(G_d,C_d);
    y = lsim(G_d*cloop_d.Si, r) + lsim(cloop_d.So*H_d, d);
    u = r - lsim(C_d, y);
    data = struct();
    data.u = u;
    data.r = r;
    data.y = y;
    data.Ts = Ts;
    data_set{itr} = data;
    parfor_progress();
end
parfor_progress(0);
%%
[mag_ori,~,wout_ori] = bode(G);
wrange = {wout_ori(1),wout_ori(end)};
mag_result_set = cell(1,max_itr);
wout_result_set = cell(1,max_itr);
parfor_progress(max_itr);
parfor itr = 1 : max_itr
    G2 = CLRIVC(data_set{itr},C_d,[na,nb,nc,nd],1);
    [mag,~,wout] = bode(G2,wrange);
    mag_result_set{itr} = mag;
    wout_result_set{itr} = wout;
    parfor_progress();
end
parfor_progress(0);

figure
for itr = 1 : max_itr
    semilogx(wout_result_set{itr}, mag2db(squeeze(mag_result_set{itr})),'r');
    hold on, box on, grid on
    drawnow
end
semilogx(wout_ori, mag2db(squeeze(mag_ori)),'b:','LineWidth',3.0);
