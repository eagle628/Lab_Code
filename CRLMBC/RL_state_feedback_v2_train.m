classdef RL_state_feedback_v2_train < handle
% CRLMBC_test
% 近似モデルとプラントが同一次元でしかできない．さらに，オブザーバんーなし
    
    properties(Constant)
        Q = diag([10,1])
        R = 1
        lambda_theta = 0.99
        lambda_omega = 0.99
        alpha = 0.001
        beta = 0.001
        gamma = 0.9
        gamma2 = 0.9
        max_episode = 5e3
    end
    
    properties
        model
        sim_N
        t
        basis_N
        c
        sigma2_basis
        sigma2_pi
    end
    
    methods
        function obj = RL_state_feedback_v2_train(model, Te, basis_N, seed)
            if nargin < 4
                seed = rng();
            end
            rng(seed)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            obj.basis_N = basis_N;
            obj.c = randn(obj.basis_N, obj.model.apx_nx);
%             obj.c = zeros(obj.basis_N, obj.model.apx_nx);
%             m = (-5:5)';
%             n = (-5:5)';
%             obj.c = [kron(m,ones(11,1)),repmat(n,11,1)]; 
            obj.sigma2_basis = 0.1*ones(obj.basis_N, 1);
%             obj.sigma2_basis = randn(obj.basis_N, 1);
            obj.sigma2_pi = 0.1;
        end
        
        function J = cost(obj, x_all, u_all)
           J = 0;
           for itr =  1 : obj.sim_N-1
               J = J + x_all(itr, :)*obj.Q*x_all(itr, :)' + u_all(itr, :)*obj.R*u_all(itr, :)';
           end
        end
        
        function R = reward(obj, x, u)
            R = -1/10*(x*obj.Q*x' + u*obj.R*u');
%             R = -1/100*(x(1)*10*x(1)' + u*obj.R*u');
        end
        
        function phi = state_basis_func(obj, x)
           phi = exp(-sum((x-obj.c).^2, 2)./(2*obj.sigma2_basis));
        end
        
        function [apx_x_all, mpc_u_all, rl_u_all, omega, theta] = actor_critic_with_eligibility_traces_episodic(obj, ini)
            rng('shuffle')
            cost_history = zeros(obj.max_episode, 1);
            reward_history = zeros(obj.max_episode, 1);
            apx_x_all = zeros(obj.sim_N,  obj.model.apx_nx);
            rl_u_all = zeros(obj.sim_N, obj.model.nu);
            mpc_u_all = zeros(obj.sim_N, obj.model.nu);
            % calculate MBC gain
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            % set episode initial
            apx_x_all(1, :) = ini';
            % params initialize
            theta = zeros(obj.basis_N, 1);
%             theta = ones(obj.basis_N, 1)*0.01;
            omega = zeros(obj.basis_N, obj.model.nu);
            for episode = 1 : obj.max_episode
                % episode initialize
                z_theta = zeros(obj.basis_N, 1);
                z_omega = zeros(obj.basis_N, 1);
                zeta = 1;
                reward = 0;
                % explration gain
                gain = 2^(-floor(episode/1000));
                for k = 1 : obj.sim_N-1
                   % MBC input
                   mpc_u_all(k, :) = -K*apx_x_all(k, :)';
                   % reinforcement learning
                   if k ~= 1
%                        r = obj.reward(apx_x_all(k, :), rl_u_all(k-1, :)+mpc_u_all(k-1, :));
%                        r = obj.reward((obj.model.apx_dynamics(apx_x_all(k, :), mu_rl+mpc_u_all(k, :)))', mu_rl+mpc_u_all(k, :));
%                        reward =  reward + obj.gamma^(k-1)*r;
                       V_k0 = obj.state_basis_func(apx_x_all(k, :))'*theta;
                       V_k1 = obj.state_basis_func(apx_x_all(k-1, :))'*theta;
                       delta = r + obj.gamma*V_k0 -V_k1;
                       z_theta = obj.gamma*obj.lambda_theta*z_theta + zeta*obj.state_basis_func(apx_x_all(k-1, :));
                       e_k1 = ((rl_u_all(k-1, :) - mu_rl)./obj.sigma2_pi)*obj.state_basis_func(apx_x_all(k-1, :));
                       z_omega = obj.gamma*obj.lambda_omega*z_omega + zeta*e_k1;
                       theta = theta + obj.alpha*delta*z_theta;
                       omega = omega + obj.beta*delta*z_omega;
                       zeta = obj.gamma2*zeta;
                   end
                   mu_rl = obj.state_basis_func(apx_x_all(k, :))'*omega;
%                    rl_u_all(k, :) = mu_rl + double(rand(1)<obj.epsilon)*obj.sigma2_pi*randn(1, obj.model.nu);
                   rl_u_all(k, :) = mu_rl + gain*obj.sigma2_pi*randn(1, obj.model.nu);
                   [apx_x_all(k+1, :), y] = obj.model.dynamics(apx_x_all(k, :), rl_u_all(k, :) + mpc_u_all(k, :));
                   r = obj.reward(apx_x_all(k+1, :), rl_u_all(k, :)+mpc_u_all(k, :));
                   reward =  reward + obj.gamma^(k-1)*r;
               end
               reward_history(episode) = reward;
               cost_history(episode) = obj.cost(apx_x_all, rl_u_all+mpc_u_all);
               callback_RL(episode, obj.t, apx_x_all, cost_history, reward_history)
            end
        end
        
        
        
        function x_all = sim(obj, ini, omega)
            x_all = zeros(obj.sim_N, obj.model.apx_nx);
            x_all(1, :) = ini';
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            for itr = 1 : obj.sim_N-1
                u_rl = obj.state_basis_func(x_all(itr, :))'*omega ;%+ obj.sigma2_pi*randn(1, obj.model.nu);
                u_mbc = -K*x_all(itr, :)';
                ne_x = obj.model.dynamics(x_all(itr,:), u_mbc+u_rl);
                x_all(itr+1, :) = ne_x';
            end
        end
        
        function [y_all] = sim_lqrcontroller(obj, ini)
            y_all = zeros(obj.sim_N, 2);
            y_all(1, :) = ini';
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            for itr = 1 : obj.sim_N-1
                y_all(itr+1, :) = obj.model.dynamics(y_all(itr,:), -K*y_all(itr,:)');
            end 
        end
    end
end
