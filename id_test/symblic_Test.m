clear all

s = sym('s');
eye_n = sym('eye_n');
K1_t = sym('K1_t');
A1 = sym('A1');
B1 = sym('B1');
L1 = sym('L1');
G2 = sym('G2');
L2 = sym('L2');
G1 = sym('G1');
eye_m = sym('eye_m');
A2 = sym('A2');
L1_t = sym('L1_t');
A1_t = sym('A1_t');


P_A = [s*eye_n-A1-K1_t*B1,L1*G2,-K1_t*B1;L2*G1,s*eye_m-A2,0;0,L1_t*G2,s*eye_m-A1_t;];
inv(P_A);

P_B =  [0;1;0;];
P_C = eye(3);

G = P_C*inv(P_A)*P_B;

n=6;m=0;
a = sym('a',[n,1]);
b = [sym('b0');sym('b',[n-m,1])];

den = (s.^(0:1:n))*[1;a];
num = (s.^(0:1:n-m))*b;
G_sym = num/den;

dG_sym = hessian(G_sym,[a;b]);
dG_sym = expand(dG_sym);
dG_sym = collect(dG_sym,'s');
params = length([a;b]);
G_func = cell(params,params,2);
for i = 1:params
    for j = 1: params
         [num_sym,den_sym] = numden(dG_sym(i,j));
          num_coe = coeffs(num_sym,'s');
          den_coe = coeffs(den_sym,'s');
          G_func(i,j,1) = {symfun(num_coe,[a;b])};
          G_func(i,j,2) = {symfun(den_coe,[a;b])};
    end
end