function FI = FI_cal(idsys, data, dim)
    N = length(data.u);
    Ts = data.Ts;

    [num, den] = tfdata(idsys.G,'v');

    tar = zeros(N, 2*dim+1);
    filter = idsys.H;

    for itr = 0 : dim-1
        num = zeros(size(num));
        num(end-itr) = 1;
        tmp = minreal(inv(filter)*sp_c2d(tf(num,den)*idsys.G, Ts), [] ,false);
        tar(:, itr+1) = lsim(tmp, data.u);
    end

    for itr = 0 : dim
        num = zeros(size(num));
        num(end-itr) = 1;
        tmp = minreal(inv(filter)*sp_c2d(tf(num,den), Ts), [], false);
        tar(:, dim+itr+1) = -lsim(tmp, data.u);
    end

    FI = tar'*tar;
end

function sys_d = sp_c2d(sys, Ts, method)
        if nargin < 3
            method = 'zoh';
        end
        [A, B, C, D] = ssdata(sys);
        nx = size(A, 1);
        nu = size(B, 2);
        no = size(C, 1);
        switch method
            case 'zoh'
                Md = expm([A, B;sparse(nu, nx+nu)]*Ts);
                Ad = Md(1:nx, 1:nx);
                Bd = Md(1:nx, nx+(1:nu));
                Cd = C;
                Dd = D;
        end
        sys_d = ss(Ad, Bd, Cd, Dd, Ts);
    end
