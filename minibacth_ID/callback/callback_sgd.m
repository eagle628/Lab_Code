function callback_sgd(x, v,  itr, y, yhat)

fprintf('%d, %e, %e\t', itr, v, x);
fprintf('\n');
opt.iteration = itr;
subplot(3,1,3), plot(itr, v, 'r*', 'MarkerSize', 5), hold on
f = gcf;
if numel(f.CurrentAxes.Children) > 100
    delete(f.CurrentAxes.Children(end));
end

callback_fmin(x, opt, [], y, yhat, []);
end

