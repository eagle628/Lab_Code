% N4SID algorithm
clear all
load('iddata_ideal')
u = w0;
y = v0;
%% constant
% system order
n = 6;
% ss estimate constant
i = 12;
h = 12;
k = 6;
j = 6;

%% Hankel Matrix
Up = hankel(u(1:i),u(i:i+j-1));
Uf = hankel(u(i+1:i+h),u(i+h:i+h+j-1));
U = [Up;Uf;];
Up_p = U(1:i+1,:);
Uf_m = U(i+2:end,:);

Yp = hankel(y(1:i),y(i:i+j-1));
Yf = hankel(y(i+1:i+h),y(i+h:i+h+j-1));
Y = [Yp;Yf;];
Yp_p = Y(1:i+1,:);
Yf_m = Y(i+2:end,:);

Wp = [Up;Yp;];
Wp_p = [Up_p;Yp_p;];

%% Oblique projections
[L,Q] = lq_de([Uf;Wp;Yf]);

L1 = [1,2];
L2 = [3,4];
L3 = [5,6];

L11 = L(1:12,L1);
L21 = L(13:36,L1);
L31 = L(37:48,L1);
L22 = L(13:36,L2);
L32 = L(37:48,L2);
L33 = L(37:48,L3);

Q1 = Q(L1,:);
Q2 = Q(L2,:);
Q3 = Q(L3,:);

%% local function

%{
L = triu(qr([U; Y]'))';
Lf = L( (2*m+l)*s+1: 2*(m+l)*s, :);
Lp = [L(1:m*s,:); L(2*m*s+1:(2*m+l)*s,:)];
Lu = L(m*s+1: 2*m*s, 1: 2*m*s);
Lfp = [Lf(:, 1: 2*m*s) - ...
(Lf(:, 1: 2*m*s)/Lu)*Lu, ...
Lf( : , 2*m*s+1: 2*(m+l)*s)];
Lpp = [Lp(:, 1: 2*m*s) - ...
(Lp(:, 1: 2*m*s)/Lu)*Lu, ...
Lp(:, 2*m*s+1: 2*(m+l)*s)];
if (norm(Lpp(: ,(2*m+l)*s-2*1 : (2*rn+l)*s), ...
'fro')) < 1e-lO
o = (Lfp*pinv(Lpp')')*Lp;
else
o = (Lfp/Lpp)*Lp;
end 
%}