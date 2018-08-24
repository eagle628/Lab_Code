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