%% NOTE
% Subspace MOESP Contiuous time Identification.
% Reference
% Identification of Continuous-time Models from Sampled Data
% y : output data
% u : output data
% d : Horizontal width (more than identificaiton dimesion)
% n : Identification dimension
% tau : prefilter parameter
%% main
function sys = moesp_CT(data, n, d, tau)
    y = data.y;
    u = data.u;
    [ndat,ny]=size(y);
    [mdat,nu]=size(u);
    if ndat ~= mdat
        error('Y and U have different length.')
    end
    if nargin < 4
        tau = data.Ts;
    elseif nargin < 3
        d = 2*n+1;
    end
    % block Hankel matrix
    N = data.N;
    Y = zeros(d*ny,N);
    U = zeros(d*nu,N);
    
    lambda = tf(1, [tau,1]);

    Y(end-ny+1:end, :) = y';
    U(end-nu+1:end, :) = u';
    t = (0:N-1)'*data.Ts;
    for s = 1 : d-1
        Y(end-(s+1)*ny+1:end-s*ny, :) = (lsim(eye(ny)*(lambda^s), y, t, [], 'foh'))';
        U(end-(s+1)*nu+1:end-s*nu, :) = (lsim(eye(nu)*(lambda^s), u, t, [], 'foh'))';
    end
    
    % LQ decomposition
    R=triu(qr([U;Y]'))';
    R=R(1:d*(ny+nu), :);

    % SVD
    R22 = R(d*nu+1:end, d*nu+1:end);
    [U1, S1, ~]=svd(R22);

    % sigular value
    singular = diag(S1);

    % C and A
    Ok = U1(:, 1:n)*diag(sqrt(singular(1:n)));
    C = Ok(1:ny, :);
    A_l = Ok(1:ny*(d-1), :)\Ok(ny+1:d*ny, :);
    A = (A_l - eye(size(A_l)))./tau;
    
    % B and D
    L1 = U1(:, n+1:end)';
    R11 = R(1:d*nu, 1:d*nu);
    R21 = R(d*nu+1:end, 1:d*nu);
    M1 = L1*R21/R11;
    m = ny*d-n;
    M = zeros(m*d, nu);
    L = zeros(m*d, ny+n);
    for k=1:d
        M((k-1)*m+1:k*m, :)=M1(:, (k-1)*nu+1:k*nu);
        L((k-1)*m+1:k*m, :)=[L1(:, (k-1)*ny+1:k*ny) L1(:, k*ny+1:end)*Ok(1:end-k*ny,:)];
    end
    DB = L\M;
    D = DB(1:ny, :);
    B = DB(ny+1:end, :)./tau;

    sys = ss(A, B, C, D);
end
