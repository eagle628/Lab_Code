function [gain_func, ang_func] = fbode(sys)
    if sys.Ts ~= 0
        sys = d2c(sys,'foh');
    end
    [num_set, den_set] = tfdata(sys);
    [no, ni] = size(num_set);

    j = sqrt(-1);
    gain_func = cell(no, ni);
    ang_func  = cell(no, ni);
    for itr1 = 1 : no
        for itr2 = 1 : ni
            num = num_set{itr1, itr2};
            den = den_set{itr1, itr2};
            G = @(w) polyval(num, j*w)/polyval(den, j*w);
            gain_func{itr1, itr2} = @(w) abs(G(w));
            ang_func{itr1, itr2}  = @(w) angle(G(w));
        end
    end
end
