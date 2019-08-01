classdef general_actor_critic_with_eligibility_traces_episodic < RL_train
    % 一般的な制御問題に対するRL_actor_criticのclass
    
    properties(Constant)
        Q = diag([10,1])
        R = 1
        lambda_theta = 0.99
        lambda_omega = 0.99
        alpha = 0.0005
        beta_mu = 0.0001
        beta_sigma = 0.01
        gamma = 0.99
        gamma2 = 0.9
        max_episode = 3.5e3
    end
    
    properties
        policy
        value
        sim_N
        t
    end
    
    methods
        function obj = general_actor_critic_with_eligibility_traces_episodic(model, Te, policy, value)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            obj.policy = policy;
            obj.value  =  value;
        end
        
        function [x_all, mpc_u_all, rl_u_all, theta_mu, w] = train(obj, ini, seed, varargin)
            if nargin < 3 || isempty(seed)
                seed = rng();
            end
            rng(seed)
            rng('shuffle')
            cost_history = zeros(obj.max_episode, 1);
            reward_history = zeros(obj.max_episode, 1);
            x_all = zeros(obj.sim_N, obj.model.true_nx);
            y_all = zeros(obj.sim_N, obj.model.ny);
            rl_u_all = zeros(obj.sim_N, obj.model.nu);
            mpc_u_all = zeros(obj.sim_N, obj.model.nu);
            % calculate MBC gain
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            tmp = strcmp(varargin, 'parallel');
            if sum(tmp)
                if  strcmp(varargin{find(tmp)+1},'off')
                    K = zeros(size(K));
                end
            end
            % set episode initial
            x_all(1, :) = ini';% When network model, local system state
            % params initialize
            % as much as possible (only information that we have)
%             w = set_w(obj, K);
            w = obj.value.get_params();
            theta_mu = obj.policy.get_params();
            theta_sigma = obj.policy.get_policy_sigma();
            % start episode learning
            for episode = 1 : obj.max_episode
                % episode initialize
                z_w = zeros(obj.value.apx_function.N, 1);
                z_theta_mu = zeros(obj.policy.apx_function.N, 1);
                z_theta_sigma = zeros(obj.model.nu, 1);
                zeta = 1;
                reward = 0;
                for k = 1 : obj.sim_N-1
                    % MBC input
                    mpc_u_all(k, :) = -K*x_all(k, :)';
%                    mpc_u_all(k, :) = 0;
                    % RL input
                    rl_u_all(k, :) = obj.policy.stocastic_policy(x_all(k, :), theta_mu);
                    %  observe S_(k+1)
                   [ne_x, y] = obj.model.dynamics(x_all(k, :)', rl_u_all(k, :) + mpc_u_all(k, :));
                    x_all(k+1, :) = ne_x';
                    y_all(k, :) = y';
                    % Get Reward r
                    if abs(x_all(k,1)) > 0.5 || abs(x_all(k,2)) > 4
                        r = -10;
                    else
                        r = obj.reward(x_all(k+1, :), rl_u_all(k, :)+mpc_u_all(k, :));
                    end
                    reward =  reward + obj.gamma^(k-1)*r;
                    % TD Erorr
                    V_k1 = obj.value.est_value(x_all(k+1, :), w);
                    V_k0 = obj.value.est_value(x_all(k, :), w);
                    delta = r + obj.gamma*V_k1 - V_k0;
                    % eligibility traces update
                    z_w = obj.gamma*obj.lambda_theta*z_w + zeta*obj.value.value_grad(x_all(k, :));
                    e_k1_mu = obj.policy.policy_grad_mu(rl_u_all(k, :), x_all(k, :), theta_mu);
                    z_theta_mu = obj.gamma*obj.lambda_omega*z_theta_mu + zeta*e_k1_mu;
                    e_k1_sigma = obj.policy.policy_grad_sigma(rl_u_all(k, :), x_all(k, :), theta_mu);
                    z_theta_sigma = obj.gamma*obj.lambda_omega*z_theta_sigma + zeta*e_k1_sigma;
                    % apx function update
                    w = w + obj.alpha*delta*z_w;
                    theta_mu = theta_mu + obj.beta_mu*delta*z_theta_mu;
                    theta_sigma = theta_sigma + obj.beta_sigma*delta*z_theta_sigma;
                    obj.policy.set_policy_sigma(theta_sigma);
                    zeta = obj.gamma2*zeta;
%                     if abs(x_all(k,1) ) > 0.5 || abs(x_all(k,2)) > 4
%                          break;
%                     end
%                     figure(3)
%                     stem(w)
%                     drawnow
                end
%                close 3
               reward_history(episode) = reward;
               cost_history(episode) = obj.cost(x_all, rl_u_all+mpc_u_all);
               figure(1)
               callback_RL(episode, obj.t, x_all, cost_history, reward_history)
               figure(2)
               [X,Y] = meshgrid(-0.5:.1:0.5, -2:.4:2);
               mesh_size = size(X, 1);
               Z = zeros(mesh_size, mesh_size);
               for itr1 = 1 : mesh_size
                   for itr2 = 1 :mesh_size
                        Z(itr1, itr2) = obj.value.est_value([X(itr1,itr2),Y(itr1,itr2)], w);
%                         Z(itr1, itr2) = -K*[X(itr1,itr2),Y(itr1,itr2)]'+obj.state_basis_func([X(itr1,itr2),Y(itr1,itr2)])'*theta_mu;
                   end
               end
               mesh(X,Y,Z)
            end
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
        
        function x_all = sim(obj, ini, w, varargin)
            x_all = zeros(obj.sim_N, obj.model.true_nx);
            x_all(1, :) = ini';
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            tmp = strcmp(varargin, 'parallel');
            if sum(tmp)
                if  strcmp(varargin{find(tmp)+1},'off')
                    K = zeros(size(K));
                end
            end
            for itr = 1 : obj.sim_N-1
%                 u_rl = obj.policy.determistic_policy(x_all(itr, :), w);
                u_rl = obj.policy.stocastic_policy(x_all(itr, :), w);
                u_mbc = -K*x_all(itr, :)';
                ne_x = (obj.model.dynamics(x_all(itr,:)', u_mbc+u_rl))';
                x_all(itr+1, :) = ne_x;
            end
        end
        
        function [x_all] = sim_lqrcontroller(obj, ini)
            x_all = zeros(obj.sim_N, 2);
            x_all(1, :) = ini';
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            for itr = 1 : obj.sim_N-1
                u_mbc = -K*x_all(itr, :)';
                x_all(itr+1, :) = (obj.model.dynamics(x_all(itr,:)', u_mbc))';
            end 
        end
        
        function w = set_w(obj, K, w0)
            if nargin < 3
                w0 = randn(obj.value.apx_function.N, 1);
            end
            [X,Y] = meshgrid(-2:.1:2, -2:.1:2);
            mesh_size = size(X, 1);
            Z = zeros(mesh_size, mesh_size);
            for itr1 = 1 : mesh_size
               for itr2 = 1 : mesh_size
                    Z(itr1, itr2) = obj.cumulative_reward(K, [X(itr1,itr2);Y(itr1,itr2)]);
               end
            end
            state = [vec(X), vec(Y)];
%             Z = -1*ones(size(Z));
            target = vec(Z);

            options = optimoptions('lsqnonlin','Display','iter','SpecifyObjectiveGradient',true);
            w = lsqnonlin(@(w)obj.apx_cost_function(state, target, w),w0,[],[],options);
% %             
% %             [X,Y] = meshgrid(-0.5:.1:0.5, -2:.4:2);
% %             mesh_size = size(X, 1);
% %             Z = zeros(mesh_size, mesh_size);
% %             for itr1 = 1 : mesh_size
% %                for itr2 = 1 :mesh_size
% %                     Z(itr1, itr2) = obj.value.est_value([X(itr1,itr2),Y(itr1,itr2)], w);
% %                end
% %             end
% %             mesh(X,Y,Z)
            
        end
        
        function [r ,apx_all, u_all] = cumulative_reward(obj, K, ini)
            apx_all = zeros(obj.sim_N, obj.model.apx_nx);
            y_all = zeros(obj.sim_N, obj.model.ny);
            u_all = zeros(obj.sim_N, obj.model.nu);
            r = 0;
            % set initial
            apx_all(1, :) = ini';
            for k = 1 : obj.sim_N-1
                u_all(k, :) = (-K*apx_all(k, :)')';
                [ne_x, y] = obj.model.approximate_dynamics(apx_all(k, :)', u_all(k, :)');
                apx_all(k+1, :) = ne_x';
                y_all(k, :) = y';
                r = r + obj.gamma^(k-1)*obj.reward(apx_all(k+1, :), u_all(k, :));
            end
        end
        
        function [error, grad] = apx_cost_function(obj, x, target, w)
            apx = zeros(size(target));
            for itr = 1 : length(target)
                apx(itr) = obj.value.est_value(x(itr, :), w);
            end
            error = apx - target;
            if nargout > 1
                grad = zeros(length(target), obj.value.apx_function.N);
                for itr = 1 : length(target)
                    grad(itr, :) = obj.value.value_grad(x(itr, :));
                end
            end
        end
    end
end

