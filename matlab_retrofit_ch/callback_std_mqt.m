function [ output_args ] = callback_std_mqt(itr,model,vals,yhat,S,lambda)
%CALLBACK_STD ?�?�?�̊֐�?�̊T?�v?�?�?�?�?�?�?�ɋL?�q
%   ?�W?�?�?�I?�ȃR?�[?�?�?�o?�b?�N?�֐�
if isfield(model,'H')&&~isempty(model.H)
    H = model.H;
else
    H = tf(1);
    H = c2d(H,1);
end
y = model.y;

vars = model.params;

for i = 1:numel(vars)
    s = vars{i};
    if startsWith(s, 'log_')
       vars{i} = s(5:end);
       vals(i) = exp(vals(i));
    end
end

if itr == 0
    fprintf('itr\t');
    %     keyboard
        fprintf('\t%-15s','f_cost');
    fprintf('\tlambda');
    for n=1:numel(model.params)
        fprintf('\t%-15s',char(vars{n}));
    end
    fprintf('\n');
    fprintf('%d\t',itr);
    fprintf('\t%.8e',S);
    fprintf('\t%.1e',lambda);
    fprintf('\t[');
    fprintf('%e\t',vals);
    fprintf(']''');
    fprintf('\n');
    if isfield(model,'H')&&~isempty(model.H)
        yf = lsim(1/H,y);
        subplot(2,1,1),plot([yf,yhat]);
        legend 'y' 'yhat'
        subplot(2,1,2),plot(yf-yhat);
    else
        subplot(2,1,1),plot([y,yhat]);
        legend 'y' 'yhat'
        subplot(2,1,2),plot([y-yhat]);
        
    end
elseif itr == inf
    fprintf('-------------------------------------------\n');
else
    fprintf('%d\t',itr);
    fprintf('\t%.8e',S);
    fprintf('\t%.1e',lambda);
    fprintf('\t[');
    fprintf('%e\t',vals);
    fprintf(']''');
    fprintf('\n');
    if isfield(model,'H')&&~isempty(model.H)
        yf = lsim(1/H,y);
        subplot(2,1,1),plot([yf,yhat], '.-');
        legend 'y' 'yhat'
        subplot(2,1,2),plot(yf-yhat, '.-');
    else
        subplot(2,1,1),plot([y,yhat], '.-');
        legend 'y' 'yhat'
        subplot(2,1,2),plot([y-yhat], '.-');
        
    end
    
    drawnow;
end


end

