function [ out ] = params2struct( vars,vals )
%PARAMS2STRUCT 推定結果を構造体に格納
%   S = params2struct(vars,vals) 

for k=1:numel(vars)
   str = char(vars{k});
   out.(str) = vals(k);
end
end

