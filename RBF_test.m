clear 
close all

basis_N = 3;
range = [-2,2];
width = (range(2)-range(1))/(basis_N-1);
m = range(1):width:range(2);
nnn = 2;
mu = m;
for itr = 1 : nnn-1
    mu = combvec(mu, m); % col vector combinater function
end
mu = mu';
sigma = 10*ones(size(mu, 1), 1);
RBF = Radial_Basis_Function(size(mu, 1), mu, sigma);

rng(6)
w = 1e-3*rand(RBF.N, 1);

func = @(x)RBF.basis_func(x);
[X,Y] = meshgrid(-4:.1:4, -4:.1:4);
mesh_size = size(X, 1);
Z = zeros(mesh_size, mesh_size);
parfor itr1 = 1 : mesh_size
   for itr2 = 1 :mesh_size
        Z(itr1, itr2) = func([X(itr1,itr2),Y(itr1,itr2)])'*w;
   end
end
mesh(X,Y,Z)
