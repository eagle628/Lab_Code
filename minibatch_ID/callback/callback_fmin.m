function out = callback_fmin(x, opt, state, y, yhat, model)

itr = opt.iteration;

if mod(itr, 1) == 0

% S = opt.fval;
% lambda = opt.lambda;
% vals = x;
if itr == 0
%     fprintf('itr\t');
%     %     keyboard
%     fprintf('\t%-15s','f_cost');
% %     fprintf('\tlambda');
%     for n=1:numel(model.params)
%         fprintf('\t%-15s',char(model.params{n}));
%     end
%     fprintf('\n');
%     fprintf('%d\t',itr);
%     fprintf('\t%.8e',S);
% %     fprintf('\t%.1e',lambda);
%     fprintf('\t[');
%     fprintf('%e\t',x);
%     fprintf(']''');
%     fprintf('\n');
    yh= yhat(x);
    subplot(3,1,1), plot([y, yh]);
    legend 'y' 'yhat'
    subplot(3,1,2), plot([y-yh]);
    drawnow
    
elseif itr == inf
%     fprintf('-------------------------------------------\n');
else
%     fprintf('%d\t',itr);
%     fprintf('\t%.8e',S);
% %     fprintf('\t%.1e',lambda);
%     fprintf('\t[');
%     fprintf('%e\t',vals);
%     fprintf(']''');
%     fprintf('\n');
    yh= yhat(x);
    subplot(3,1,1), plot([y, yh]);
    legend 'y' 'yhat'
    subplot(3,1,2), plot([y-yh]);
    
    drawnow;
end
end
out = false;
end

