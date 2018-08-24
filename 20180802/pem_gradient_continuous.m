function [dyhat,sys_update] = pem_gradient_continuous(sys, u, t)
    num = sys.Numerator{:};
    den = sys.Denominator{:};
    params = numel(num)+numel(den)-1;
    dim = numel(den)-1;
    dyhat = zeros(length(u),params);
    for itr = 1 : numel(num)
        num = num*0;
        num(itr) = 1;
        dsys = tf(num,den);
        dyhat(:,dim+itr) = lsim( dsys, u, t);
    end
    for itr = 1 : dim
        num = num*0;
        num(itr+1) = 1;
        dsys = -tf(num,den)*sys;
        dyhat(:,itr) = lsim( dsys, u, t);
    end
    den(2:end) = den(2:end) - mean(dyhat(:,1:dim));
    num = sys.Numerator{:} - mean(dyhat(:,dim+1:end));
    sys_update = tf(num,den);
end
