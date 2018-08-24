clear
%% Load Data
% y1 Add node
% y2 Add controller
% y3 origianl
% y4 Remove controller
load('Response.mat')
Ts = mean(diff(t_s));
ch_t = ch_t/Ts+1;
Fs = 1/Ts;
%% FFT
state = 2;
[ P1_1, P2_1, L] = fft_func( y1(ch_t:length(y1),2), 3);
[ P1_2, P2_2] = fft_func( y2(ch_t:length(y2),2), 3);
[ P1_3, P2_3] = fft_func( y3(ch_t:length(y3),2), 3);
[ P1_4, P2_4] = fft_func( y4(ch_t:length(y4),2), 3);

figure
f = Fs*(0:(L/2))/L;
Peak = max(P1_3);
plot(f,P1_3,'g-','linewidth',3.0)
hold on
grid on
plot(f,P1_1,'b-','linewidth',0.8)
plot(f,P1_2,'r-','linewidth',0.8)
plot(f,P1_4,'k-','linewidth',0.8)
title('Single-Sided Amplitude Spectrum')
xlabel('f (Hz)')
xlim([0 1])
legend('Original(Extend2)','AddNode(Extend2)','Original(Extend3)','Original(Extend1)')

%% local function
function [ P1, P2, L] = fft_func( signal, window)
    L = length(signal);
    %1.矩形窓 2.Blackman Window 3.Hamming Window  4.Hann Window
    switch window
        case 1
            w = 1;
        case 2    
            w = blackman(L);
        case 3   
            w = hamming(L);
        case 4  
            w = hann(L);
    end
    Y = fft(signal.*window);
    P2 = abs(Y/L);    % Both side amplitude
    P1 = P2(1:L/2+1); % single side amplitude
    P1(2:end-1) = 2*P1(2:end-1);
end
