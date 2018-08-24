%{
coeffs は次数が低い順に係数を取り出す．つまり，反転しないといけない． 
%}
close all
clear

Ts = 0.1;
dim = 6;

P = rss(dim);
P = tf(P);
P = minreal(P);

a = sym('a',[dim, 1]);
% assume(a, 'real')
b = sym('b',[dim+1, 1]);
% assume(b, 'real')

syms s

Phat = s.^(dim:-1:0)*b/(s.^(dim:-1:0)*[1;a]);
dPhat = cell( size(a,1)+size(b,1), 1);
for itr1 = 1 : size(a,1)
    dPhat{itr1} = diff(Phat,a(itr1));
end
for itr2 = 1 : size(b,1)
    dPhat{itr1+itr2} = diff(Phat,b(itr2));
end    



% [numP, denP] = tfdata(P, 'v');
% P2 = subs(Phat, [a ; b], [denP(2:end), numP]');
% [num_s, den_s] = numden(P2);
% num = double(fliplr(coeffs(num_s, s)))
% den = double(fliplr(coeffs(den_s, s)))
% KKK = tf(num, den);
% KK = minreal(KKK)
