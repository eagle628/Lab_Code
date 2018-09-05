function callback_sgd(x, v,  itr, y, yhat)

fprintf('%d, %e, %e\t', itr, v, x);
fprintf('\n');
opt.iteration = itr;
callback_fmin(x, opt, [], y, yhat, []);
end

