%% NOTE
% Identification of Continuous-time Models from Sampled Data.pdf 
% Chapter 5
% section 5.5.3
% Two Step CLIVC2 algorithm
%% main
clear
close all
%% Prepare Plan and Controller and Noise Filter
% Plant
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
% Controller
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
% Noise filter
HHH = 1;
rng(5)
while true
    H_d = tf([1,randn(1,HHH)],[1,randn(1,HHH)],Ts);
    if isstable(H_d)
        break;
    end
end

% [~, den] = tfdata(G_d,'v');
% H = tf(1, den);
% H_d = c2d(H, Ts, 'foh');
H_d = tf(1,1,Ts);

nc = HHH;
nd = HHH;

%% Generate Response : Availabel r1 r2 y u
max_itr = 50;
data_set = cell(1,max_itr);
rng('shuffle')

for itr = 1 : max_itr
    r1 = randn(N, 1);
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
end
%% 

% % [G1,G2] = CLIVC2(data,C_d,[na,nb],1,true);
% % figure('name','Book Method')
% % bode(G,'b',G1,'r',G2,'g')
% % legend('original','First','Second')
% [G1,G2] = CLIVC2(data,C_d,[na,nb],1,false);
% figure('name','My Method')
% bode(G,'b',G1,'r',G2,'g')
% legend('original','First','Second')
[mag,~,wout] = bode(G,{1e-3,1e1});
figure
semilogx(wout, mag2db(squeeze(mag)),'b:','LineWidth',3.0);
hold on
for itr = 1 : max_itr
    G2 = CLIVC2(data_set{itr},C_d,[na,nb],0.1);
    [mag,~,wout] = bode(G2,{1e-3,1e1});
    semilogx(wout, mag2db(squeeze(mag)),'r');
    drawnow
end
