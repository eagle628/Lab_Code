%% NOTE
% % Calculating Bode fucnttion of Continuous sytem
% % gain and angle are calculating function-handle
% % So, You need to "fplot" drawing bode plot.
% % EX)
% %    [MAG,PHASE] = fbode(sys);
% %    fplot(@(w)mag2db(MAG{1}(w)),[wmin, wmax])
% %    fplot(@(w)PHASE{1}(w),[wmin, wmax])
% % Only Minimum Phas system
%% main
function [gain_func, ang_func] = fbode(sys)
    if sys.Ts ~= 0
        sys = d2c(sys,'foh');
    end
    [num_set, den_set] = tfdata(sys);
    [no, ni] = size(num_set);

    gain_func = cell(no, ni);
    ang_func  = cell(no, ni);
    for itr1 = 1 : no
        for itr2 = 1 : ni
            num = num_set{itr1, itr2};
            den = den_set{itr1, itr2};
            G = @(w) polyval(num, 1i*w)./polyval(den, 1i*w);
            gain_func{itr1, itr2} = @(w) abs(G(w));
            ang_func{itr1, itr2}  = @(w) rad2deg(angle(G(w)));
        end
    end
end
