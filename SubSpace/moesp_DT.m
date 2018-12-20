%% NOTE
% Subspace MOESP Discrete time Identification.
% Reference
% Tohru Katayama, "Subspace Methods for System Identification", Springer
% 2005.
% y : output data
% u : output data
% d : Horizontal width (more than identificaiton dimesion)
% n : Identification dimension
%% main
function sys = moesp_DT(data, n, d)
    y = data.y;
    u = data.u;
    [ndat,ny]=size(y);
    [mdat,nu]=size(u);
    if ndat ~= mdat
        error('Y and U have different length.')
    end
    if nargin < 3
        d = 2*n+1;
    end

    % block Hankel matrix
    N = ndat-d+1;
    Y = zeros(d*ny,N);
    U = zeros(d*nu,N);
%     sN = sqrt(N);
%     sy = y'/sN;
%     su = u'/sN;
    sy = y';
    su = u';
    for s=1:d
        Y((s-1)*ny+1:s*ny, :)=sy(:, s:s+N-1);
        U((s-1)*nu+1:s*nu, :)=su(:, s:s+N-1);
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
    C=Ok(1:ny, :);
    A=Ok(1:ny*(d-1), :)\Ok(ny+1:d*ny, :);

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
    DB=L\M;
    D=DB(1:ny, :);
    B=DB(ny+1:end, :);

    sys = ss(A, B, C, D, data.Ts);
end
