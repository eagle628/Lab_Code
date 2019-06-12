classdef LQR_model_free_state_feedback_train < handle
    %UNTITLED このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Constant)
        Q = diag([10,1])
        R = eye(1)
        max_iteration = 1e3
        sigma2 = 1e-4
        epsilon = 1e-4
    end
        
    properties
        model                  
        sim_N
        t
        K0
    end
   
    
    
    methods
        function obj = LQR_model_free_state_feedback_train(model, Te)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            obj.K0 = [];
        end
        
        function set_initial_K(obj, seed)
            if nargin < 2
                seed = 0;
            end
            rng(seed)
            obj.K0 = dlqr(obj.model.A, obj.model.B, eye(size(obj.model.A)), eye(size(obj.model.B,2)));
            rng('shuffle')
        end
        
        function [x_all, u_all] = sim(obj, ini, K)
            x_all = zeros(obj.sim_N, 2);
            x_all(1, :) =  ini';
            u_all = zeros(obj.sim_N, 1);
            for itr1 = 1 : obj.sim_N-1
                u_all(itr1, :) = -K*x_all(itr1, :)';
                ne_x = obj.model.dynamics(x_all(itr1, :), u_all(itr1, :));
                x_all(itr1+1, :) = ne_x';
            end
        end
        
        function  K = train(obj,ini)
            seed = 6;
            rng(seed)
            if isempty(obj.K0)
                obj.set_initial_K();
            end
            x_all = zeros(obj.sim_N, 2);
            x_all(1, :) =  ini';
            u_all = zeros(obj.sim_N, 1);
            K = obj.K0;
            for itr1 = 1 : obj.sim_N-1
                u_all(itr1, :) = -K*x_all(itr1, :)'+ obj.sigma2*randn(1);
                ne_x = obj.model.dynamics(x_all(itr1, :), u_all(itr1, :));
                x_all(itr1+1, :) = ne_x';
            end
            [delta_xx, eye_xx, eye_xu] = intermeditate_cal(x_all, u_all);
            for itr2 = 1: obj.max_iteration
                [P, K] = updator(delta_xx, eye_xx, eye_xu, K, obj.Q, obj.R);
                if ~exist('P_set')
                    P_set = zeros(size(P, 1), obj.max_iteration);
                    P_set(:, itr2) = P;
                else
                    P_set(:, itr2) = P;
                    if norm(P_set(:, itr2) - P_set(:, itr2-1)) < obj.epsilon
                        break;
                    end
                end
            end
            rng('shuffle')
        end
    end
    
end

%% local
function [P, K] = updator(delta_xx, eye_xx, eye_xu, K, Q, R)
    nx = size(delta_xx, 2);
    nu = size(eye_xu, 2)/nx;
    Theta = [delta_xx, -2*eye_xx*(kron(eye(nx), K'*R)) - 2*eye_xu*kron(eye(nx), R)];
    Sigma = -eye_xx*vec(Q);
    P_K = (Theta'*Theta)^(-1)*Theta'*Sigma;
    k_n = nu*nx;
    K = P_K(end-k_n+1:end);
    K = transpose(vec2mat(K, nu));
    P = P_K(1:end-k_n+1);
end

function [delta_xx, eye_xx, eye_xu] = intermeditate_cal(x, u)
    delta_xx = diff(x, 1, 1);
    [mx, nx] = size(x);
    kron_xx_tmp = zeros(mx, nx*nx);
    for itr = 1 : nx*nx
        kron_xx_tmp(:, itr) = x(:, ceil(itr/nx)).*x(:, mod(itr-1, nx)+1);
    end
%                     kron_xx_tmp = kron_xx_tmp * 0.01;
    eye_xx = kron_xx_tmp(1:end-1, :) + kron_xx_tmp(2:end, :);
    nu = size(u, 2);
    kron_xu_tmp = zeros(mx, nx*nu);
    for itr = 1 : nx*nu
        kron_xu_tmp(:, itr) = x(:, ceil(itr/mx)).*u(:, mod(itr-1, nu)+1);
    end
%                     kron_xu_tmp = kron_xu_tmp*0.01;
    eye_xu = kron_xu_tmp(1:end-1, :) + kron_xu_tmp(2:end, :);
end

% 数値積分は適当な台形則での近似
     
