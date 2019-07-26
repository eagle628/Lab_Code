close all
clear

model = CRLMBC_test_model(0.5,0.15,9.8,0.5,0.1);

basis_N = 21^2;

nnn = sqrt(basis_N);
range = [-5,5];
width = (range(2)-range(1))/(nnn-1);
m = (range(1):width:range(2))';
c = [kron(m,ones(nnn,1)),repmat(m,nnn,1)]; 
sigma2_basis = 0.25*ones(basis_N, 1);

% cumulutive 
Q = diag([10,1]);
R = diag(1);
K = dlqr(model.A, model.B, Q, R);
gamma = 0.9;

Te = 3;
ini = [1;1];
r = cumulative_reward(model, K, gamma, Te, ini, Q, R);
figure('Name','True Value function')
[X,Y] = meshgrid(-1:.1:1, -1:.1:1);
mesh_size = size(X, 1);
Z = zeros(mesh_size, mesh_size);
for itr1 = 1 : mesh_size
   for itr2 = 1 :mesh_size
        Z(itr1, itr2) = cumulative_reward(model, K, gamma, Te, [X(itr1,itr2);Y(itr1,itr2)], Q, R);
   end
end
meshc(X,Y,Z)

state = [vec(X), vec(Y)];
target = vec(Z);

w0 = zeros(basis_N, 1);

options = optimoptions('lsqnonlin','Display','iter');
w = lsqnonlin(@(w)apx_cost_function(state, c, sigma2_basis, target, w),w0,[],[],options);


figure('Name','apx_function')
[X2,Y2] = meshgrid(-1:.1:1, -1:.1:1);
mesh_size = size(X2, 1);
Z2 = zeros(mesh_size, mesh_size);
for itr1 = 1 : mesh_size
   for itr2 = 1 :mesh_size
        Z2(itr1, itr2) = state_basis_func([X(itr1,itr2),Y(itr1,itr2)], c, sigma2_basis)'*w;
   end
end
meshc(X2,Y2,Z2)






%% local
function [r ,x_all, u_all] = cumulative_reward(model, K, gamma, Te, ini, Q, R)
    sim_N = Te/model.Ts + 1;
    x_all = zeros(sim_N, model.apx_nx);
    u_all = zeros(sim_N, model.nu);
    r = 0;
    % set initial
    x_all(1, :) = ini';
    for k = 1 : sim_N-1
        u_all(k, :) = (K*x_all(k, :)')';
        ne_x = (model.A-model.B*K)*x_all(k, :)';
        x_all(k+1, :) = ne_x';
        r = r + gamma^(k-1)*reward(x_all(k+1, :), u_all(k, :), Q, R);
    end
end

function R = reward(x, u, Q, R)
    R = -1/10*(x*Q*x' + u*R*u');
end

function phi = state_basis_func(x, c, sigma2_basis)
   phi = exp(-sum((x-c).^2, 2)./(2*sigma2_basis));
end

function error = apx_cost_function(x, c, sigma2_basis, target, w)
    apx = zeros(size(target));
    for itr = 1 : length(target)
        apx(itr) = state_basis_func(x(itr, :), c, sigma2_basis)'*w;
    end
    error = apx - target;
end
