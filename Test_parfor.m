a = 1:1000;
b = zeros(numel(a),1);
parfor i = 1 : numel(a)
    b(i) = a(i);
    disp(i)
end