% Householder
function a=householder(A)
    [n,n] = size(A);
    a = A;
    for k=1:n-1
        b = a(k+1:n,k);
        B = a(k+1:n,k+1:n);
        s = - sign(b(1)) * norm(b);
        w = b;
        w(1) = w(1) - s;
        normw2 = 2 * s * (s - b(1));
        p = (B * w) / (normw2/2);
        alpha = (w'*p) / normw2;
        q = p - alpha * w;
        a(k+1:n,k+1:n) = B - w * q' - q * w';
        a(k+2:n,k) = zeros(n-k-1,1);
        a(k,k+2:n) = zeros(n-k-1,1)';
        a(k+1,k) = s;
        a(k,k+1) = s;
    end
end 
