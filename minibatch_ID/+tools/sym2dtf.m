function [ tfobj ] = sym2dtf( symobj, vars, varvals, Ts)
%SYM2TF Convert symbolic math rationals to transfer function
% [ tfobj ] = sym2tf( symobj, vars, varvals)
%     this function perform subs(symobj,vars,sym(varvals,'d'))

if ~isa(symobj,'sym')
    tfobj=tf(symobj);
    return;
end

if nargin==4
    symobj=subs(symobj,vars,sym(varvals,'d'));
end
% if symobj==0
%     tfobj=tf(0);
%     return;
% end
if ~isa(symobj,'sym')
    tfobj=tf(symobj);
    return;
end


[n,d]=numden(symobj);
dc=coeffs(d);
n=n/dc(1);
d=d/dc(1);
num=sym2poly(n);
den=sym2poly(d);
fact=max(abs(den));
num=num/fact;
den=den/fact;

tfobj=tf(num,den,Ts);
end