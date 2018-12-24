function sys = CLRIVC(data, C, condition, lambda)

    % load data
    u  = data.u;
    y  = data.y;
    Ts = data.Ts;
    N  = length(u);
    % contoroller
    if C.Ts == 0
        C_d = c2d(C, Ts, 'foh');
    else
        C_d = C;
    end
    if isfield(data, 'r')
        r = data.r;
    else
        % clculate r
        r = u + lsim(C_d, y);
    end
    % conditon
    na = condition(1);
    nb = condition(2);
    nc = condition(3);
    nd = condition(4);

    % chose Initial CTfilter : f_c(p)
    if nargin < 4
        lambda = 1;
    end
    beta = 1;
    s = tf('s');
%     init_sys = (beta/(s+lambda))^na;
    init_sys = (beta/(lambda*s+1))^na;
    [num, den] = tfdata(init_sys,'v');
    % Iteration Part
    for itr1 = 1 : 3
        % step : 1
        y_f_c = zeros(N, na+1);
        for itr2 = 0 : na
            num = zeros(size(num));
            num(end-itr2) = 1;
            y_f_c(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), y);
        end
        u_f_c = zeros(N, nb);
        for itr2 = 0 : nb-1
            num = zeros(size(num));
            num(end-itr2) = 1;
            u_f_c(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), u);
        end
        phi_f_c = [-y_f_c(:,2:end), u_f_c];
        rho = (phi_f_c'*phi_f_c)\(phi_f_c'*y_f_c(:,1));
        % step : 2(a)
        G1 = tf((rho(na+1:end))',[1,(rho(1:na))']);
        G1_d = c2d(G1, Ts, 'foh');
        cloop_d_1 = loopsens(G1_d,C_d);
        % calculate noise free
        x = lsim(G1_d*cloop_d_1.Si, r);
        nu = lsim(cloop_d_1.Si, r);
        % Eq 5.39
        x_f_c = zeros(N, na);
        for itr2 = 0 : na-1
            num = zeros(size(num));
            num(end-itr2) = 1;
            x_f_c(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), x);
        end
        % Eq 5.40
        nu_f_c = zeros(N, nb);
        for itr2 = 0 : nb-1
            num = zeros(size(num));
            num(end-itr2) = 1;
            nu_f_c(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), nu);
        end
        zeta_f_c = [-x_f_c, nu_f_c];
        if itr1 == 1
            zeta_f_clrivc = zeta_f_c;
            phi_f_clrivc = phi_f_c;
            y_f_clrivc = y_f_c(:, 1);
        else
            % Step :2(b)
            xi = y - x;
            noise_sys = armax(iddata(xi,[],Ts),[nc,nd]);
            % Step :2(c)
            noise_sys = tf(noise_sys.A,noise_sys.C,Ts);
            % Step :2(d)
            zeta_f_clrivc = zeros(size(zeta_f_c));
            for itr2 = 1 : size(zeta_f_c, 2)
                zeta_f_clrivc(:, itr2) = lsim(noise_sys, zeta_f_c(:,itr2));
            end
            phi_f_clrivc = zeros(size(phi_f_c));
            for itr2 = 1 : size(phi_f_c, 2)
                phi_f_clrivc(:, itr2) = lsim(noise_sys, phi_f_c(:,itr2));
            end
            y_f_clrivc = lsim(noise_sys, y_f_c(:,1));
        end

        rho = (zeta_f_clrivc'*phi_f_clrivc)\(zeta_f_clrivc'*y_f_clrivc);
        G2 = tf((rho(na+1:end))',[1,(rho(1:na))']);
        [num, den] = tfdata(G2,'v');
    end
    sys = struct();
    sys.G = G2;
    sys.H = inv(noise_sys);
end
