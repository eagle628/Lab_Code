function [v, dv, dv_num] = test_diff(func, theta, h, pararell)
[v, dv] = func(theta);
if nargin < 3
    h = 1e-8;
end
if nargin < 4
    if islogical(h)
        pararell = h;
        h = 1e-8;
    else
        pararell  = false;
    end
end
dv_num = zeros(size(v, 1), size(v,2), numel(theta));
parfor_progress(numel(theta));
if pararell
    parfor itr = 1:numel(theta)
        theta1 = theta;
        theta1(itr) = theta1(itr) + h;
        v1 = func(theta1);
        dv_num(:, :, itr) = (v1-v)/h;
        parfor_progress();
    end
else
    for itr = 1:numel(theta)
        theta1 = theta;
        theta1(itr) = theta1(itr) + h;
        v1 = func(theta1);
        dv_num(:, :, itr) = (v1-v)/h;
        parfor_progress();
    end
end
parfor_progress(0);
dv_num = squeeze(dv_num);
end

