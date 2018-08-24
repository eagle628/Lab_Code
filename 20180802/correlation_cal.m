function [correlation,lgs] = correlation_cal(y, dy, maxlag)
    [ ~, colum] = size(dy);
%     correlation = zeros(row*2-1,colum);
    correlation = zeros(colum,maxlag*2+1);
    for itr = 1 : colum
        cr = xcorr( y, dy(:,itr), maxlag,'coeff');
        correlation(itr,:) = cr;
    end 
    if nargout > 1
        lgs = -maxlag:1:maxlag;
    end
end
