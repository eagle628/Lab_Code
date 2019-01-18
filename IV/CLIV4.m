function sys = CLIV4(data, C_d, n)
    % load data
    u  = data.u;
    y  = data.y;
    Ts = data.Ts;
    N  = length(u);
    % contoroller
    if C_d.Ts == 0
        C_d = c2d(C_d, Ts);
    end
    if isfield(data, 'r')
        r = data.r;
    else
        % clculate r
        r = u + lsim(C_d, y);
    end
    % conditon
%     n = condition(1);

    % IV method
    % first stage
    phi = zeros(N+n-1, 2*n);
    for itr = n : -1 :1
        tmp_st = n-itr;
        phi(tmp_st+1:tmp_st+N, itr) = -y;
    end
    for itr = n : -1 :1
        tmp_st = n-itr;
        phi(tmp_st+1:tmp_st+N, n+itr) = u;
    end   
    theta1 = (phi'*phi)\(phi'*[y;zeros(n-1,1)]);
    G1_d = tf(theta1(n+1:end)',theta1(1:n)', Ts, 'Variable','q^-1');
    % Second stage
    Loop = loopsens(G1_d, C_d);
    tilde_y1 = lsim(Loop.Ti, r);
    tilde_u1 = lsim(Loop.Si, r);
    z1 = zeros(N+n-1, 2*n);
    for itr = n : -1 :1
        tmp_st = n-itr;
        z1(tmp_st+1:tmp_st+N, itr) = -tilde_y1;
    end
    for itr = n : -1 :1
        tmp_st = n-itr;
        z1(tmp_st+1:tmp_st+N, n+itr) = tilde_u1;
    end
    theta2 = (z1'*phi)\(z1'*[y;zeros(n-1,1)]);
    G2_d = tf(theta2(n+1:end)',theta2(1:n)',  Ts, 'Variable','q^-1');
    % Third Stage
    w = lsim(tf(theta2(n+1:end)', 1, Ts, 'Variable','q^-1'), y) ...
            - lsim(tf(theta2(1:n)', 1, Ts, 'Variable','q^-1'), u);
    L = armax(iddata(w,[],Ts),[n,n]);
    % Fourth Stage
    Loop = loopsens(G2_d, C_d);
    tilde_y2 = lsim(Loop.Ti, r);
    tilde_u2 = lsim(Loop.Si, r);
    z2 = zeros(N+n-1, 2*n);
    for itr = n : -1 :1
        tmp_st = n-itr;
        z2(tmp_st+1:tmp_st+N, itr) = -tilde_y2;
    end
    for itr = n : -1 :1
        tmp_st = n-itr;
        z2(tmp_st+1:tmp_st+N, n+itr) = tilde_u2;
    end
    % filtering estimate noise sys
    phi_T = zeros(size(phi));
    for itr = 1 : 2*n
        phi_T(:, 1) = lsim(L, phi(:, 1));
    end
    y_T = lsim(L, y);
    theta_cliv4 =  (z2'*phi_T)\(z2'*[y_T;zeros(n-1,1)]);
    G_cliv4_d = tf(theta_cliv4(n+1:end)',theta_cliv4(1:n)', Ts, 'Variable','q^-1');
    
    sys = struct();
    sys.G = G_cliv4_d;
    sys.H = L;
end
