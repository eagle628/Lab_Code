clear 
close all
%% Set Transfer Function
s = tf('s');
z = tf('z');

Ps = (s+1)/(s+2);
%Ps = 1/(s+2);
Cs = 2*(s+4)/(s+3);
%Cs = 2/(s+3);
Hs = 1;
Vs = 1;

Ts = 0.1;
P = c2d(Ps,Ts);
C = c2d(Cs,Ts);
H = Hs;
V = Vs;

%% Bode
min = -1;
max =  1;
[mag_P ,ph_P,w_P]   = bode(P,{10^min,10^max});
[mag_Ci ,ph_Ci,w_Ci] = bode(C,{10^min,10^max});

fig1 = figure;
ax1 = subplot(2,1,1);
hold on;
ax2 = subplot(2,1,2);
hold on;

plot(ax1,w_P,mag2db(reshape(mag_P(1,1,:),[],1)),'k:','linewidth',3)
plot(ax1,w_Ci,mag2db(reshape(mag_Ci(1,1,:),[],1)),'b:','linewidth',3)


plot(ax2,w_P,reshape(ph_P(1,1,:),[],1),'k:','linewidth',3)
plot(ax2,w_Ci,reshape(ph_Ci(1,1,:),[],1),'b:','linewidth',3)

ax1.XScale = 'log';
ax1.XMinorGrid = 'on';
ax1.YMinorGrid = 'on';
ax1.Box = 'on';
ax1.XLabel.String = 'Frequency [rad/s]';
ax1.YLabel.String = 'Gain [dB]';

ax2.XScale = 'log';
ax2.XMinorGrid = 'on';
ax2.YMinorGrid = 'on';
ax2.Box = 'on';
ax2.XLabel.String = 'Frequency [rad/s]';
ax2.YLabel.String = 'Phase [deg]';

%% identification
G_eu = C*H/(1-C*P);
G_ey = H/(1-C*P);

G_ru = V/(1-C*P);
G_ry = P*V/(1-C*P);

power = [ 0, 0.1, 1.0, 10];
line = {'g-','g:','g--','g-.'};
for itr = 1:numel(power)
    N = 1000;
    rng(10);
    e = rand(N,1)*power(itr);
    rng(6)
    r = rand(N,1)*1;
    t = (0:N-1)'*Ts;
    u = lsim(G_eu,e,t)+lsim(G_ru,r,t);
    y = lsim(G_ey,e,t)+lsim(G_ry,r,t);
    data = iddata(y,u,Ts);
    nk = 0;
    model_armax = armax(data,[1,2,1,nk]);
    [mag_id ,ph_id,w_id] = bode(model_armax,{10^min,10^max});
    plot(ax1,w_id,mag2db(reshape(mag_id(1,1,:),[],1)),line{itr},'linewidth',1.5)
    plot(ax2,w_id,reshape(ph_id(1,1,:),[],1),line{itr},'linewidth',1.5)
end


legend(ax2,'G_{21}','G_{12}^{-1}','id','location','southwest')
