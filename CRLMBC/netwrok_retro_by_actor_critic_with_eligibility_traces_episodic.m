classdef netwrok_retro_by_actor_critic_with_eligibility_traces_episodic < RL_train
    % レトロフィット制御問題に対するRL_actor_criticのclass
    
    properties(Constant)
        Q1 = 1
        Q = diag([1000,1])
        R = 1
        lambda_theta = 0.99
        lambda_omega = 0.99
        alpha = 0.005
        beta = 0.001
        alpha_f = 0.001
        gamma = 0.9
        gamma2 = 0.9
        max_episode = 3.5e3
    end
    
    properties
        sim_N
        t
    end
    
    methods
        function obj = netwrok_retro_by_actor_critic_with_eligibility_traces_episodic(model, Te, basis_N, seed)
            rng('shuffle')
            if nargin < 4
                seed = rng();
            end
            rng(seed)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            range = [-1,1];
            width = (range(2)-range(1))/(basis_N-1);
            m = range(1):width:range(2);
            nnn = obj.model.local_nx + obj.model.rect_nx;
            mu = m;
            for itr = 1 : nnn-1
                mu = combvec(mu, m); % col vector combinater function
            end
            mu = mu';
            sigma = 0.5*ones(size(mu, 1), 1);
            RBF1 = Radial_Basis_Function(size(mu, 1), mu, sigma);
            sigma_pi = 1;
            obj.policy = policy_RBF(RBF1, sigma_pi);
            obj.value  =  value_RBF(RBF1);
        end
        
        function [local_x_all, mpc_u_all, rl_u_all, theta_mu, w, reward_history] = train(obj, ini)
            rng('shuffle')
            reward_history = zeros(obj.max_episode, 1);
            local_x_all = zeros(obj.sim_N, obj.model.local_nx);
            env_x_all = zeros(obj.sim_N, obj.model.env_nx);
            rect_x_all = zeros(obj.sim_N, obj.model.rect_nx);
            y_xhat_w_v_all = zeros(obj.sim_N, obj.model.ny+obj.model.ny+obj.model.nw+obj.model.nv);
            rl_u_all = zeros(obj.sim_N, obj.model.nu);
            mpc_u_all = zeros(obj.sim_N, obj.model.nu);
            % calculate MBC gain
            if obj.model.apx_nx == 0
%                 K = lqr(obj.model.A, obj.model.B, obj.Q, obj.R);
                obj.model.set_controlled_system(obj.Q, obj.R);
            else
%                 K = lqr(obj.model.A, obj.model.B, blkdiag(obj.Q, eye(obj.model.apx_nx)*obj.Q1), obj.R);
                obj.model.set_controlled_system(blkdiag(obj.Q, eye(obj.model.apx_nx)*obj.Q1), obj.R);
            end
%             K1 = K(1:obj.model.local_nx);
%             K2 = K(obj.model.local_nx+1 : end);
            % set episode initial
            local_x_all(1, :) = ini';% When network model, local system state
            % local noise
%             d_L = zeros(obj.sim_N, 2);
            d_L = randn(obj.sim_N, 2);
            % params initialize
            % as much as possible (only information that we have)
%             w = set_w(obj, K);
            w = obj.value.get_params();
            theta_mu = obj.policy.get_params();
            % start episode learning
            for episode = 1 : obj.max_episode
                % episode initialize
                z_w = zeros(obj.value.apx_function.N, 1);
                z_theta_mu = zeros(obj.policy.apx_function.N, 1);
%                 z_theta_sigma = zeros(obj.basis_N, 1);
                zeta = 1;
                reward = 0;
                % explration gain
%                 gain = 2^(-floor(episode/1000));
                gain = 1;
%                 gain = 1;
                % initialize fisher matrix
%                 G = eye(obj.policy.apx_function.N);
                for k = 1 : obj.sim_N-1
                   % MBC input
                   mpc_u_all(k, :) = 0; % LQR Controller already implemented
%                    mpc_u_all(k, :) = ([K1,-K2]*rect_x_all(k, :)')' -(K1*ywv_all(k, 1:obj.model.ny)')';
                   % RL input
                   rl_u_all(k, :) = obj.policy.stocastic_policy([local_x_all(k, :), rect_x_all(k, :)], theta_mu);
                   %  observe S_(k+1)
                   [local_ne_x, env_ne_x, ne_ywv, rect_ne_x] = ...
                        obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)',rect_x_all(k, :)', y_xhat_w_v_all(k, :)', rl_u_all(k, :) + mpc_u_all(k, :), d_L(k, :)');
                   
                    local_x_all(k+1, :) = local_ne_x';
                   env_x_all(k+1, :) = env_ne_x';
                   y_xhat_w_v_all(k+1, :) = ne_ywv';
                   rect_x_all(k+1, :) = rect_ne_x';
                   % Get Reward r
                   if abs(local_x_all(k,1) ) > 1
                       r = -100;
                   else
                       r = obj.reward(rect_x_all(k+1, :), rl_u_all(k, :)+mpc_u_all(k, :));
%                         r = 0;
                   end
                   reward =  reward + obj.gamma^(k-1)*r;
                   % TD Erorr
                   V_k1 = obj.value.est_value([local_x_all(k+1, :), rect_x_all(k+1, :)], w);
                   V_k0 = obj.value.est_value([local_x_all(k, :), rect_x_all(k, :)], w);
%                     V_k1 = 0;V_k0 = 0;
                   delta = r + obj.gamma*V_k1 - V_k0;
                   % eligibility traces update
                   z_w = obj.gamma*obj.lambda_theta*z_w + zeta *obj.value.value_grad([local_x_all(k, :), rect_x_all(k, :)]);
                   e_k1_mu = obj.policy.policy_grad(rl_u_all(k, :), [local_x_all(k, :), rect_x_all(k, :)], theta_mu);
%                    G = 1/(1-obj.aplha_f)*(G - obj.alpha_f*(G*e_k1_mu)*(G*e_k1_mu)'/(1-obj.alpha_f+obj.alpha_f*e_k1_mu'*G*e_K1_mu));
                   z_theta_mu = obj.gamma*obj.lambda_omega*z_theta_mu + zeta*e_k1_mu;
%                        e_k1_sigma = ((rl_u_all(k-1, :) - mu_rl).^2/(pi_sigma^2)-1)*obj.state_basis_func2(apx_x_all(k-1, :));
%                        z_theta_sigma = obj.gamma*obj.lambda_omega*z_theta_sigma + zeta*e_k1_sigma;
                   % apx function update
                   w = w + obj.alpha*delta*z_w;
                   theta_mu = theta_mu + obj.beta*delta*z_theta_mu;
                   zeta = obj.gamma2*zeta;
%                     if abs(apx_x_all(k,1) ) > 1
%                         break;
%                     end
               end
                reward_history(episode) = reward;
                subplot(3,1,1)
                plot(obj.t, local_x_all)
                title(['Episode-',num2str(episode)])
                ylabel('z')
                grid on
                subplot(3,1,2)
                plot(obj.t, rect_x_all)
                grid on
                ylabel('\hat{z}')
                subplot(3,1,3)
                plot(nonzeros(reward_history),'-b')
                ylabel('Culumative Reward')
                drawnow
%                cost_history(episode) = obj.cost([local_x_all(k, :), rect_x_all(k, :)], rl_u_all+mpc_u_all);
%                figure(1)
%                callback_RL(episode, obj.t, [local_x_all(k, :), rect_x_all(k, :)], cost_history, reward_history)
            end
        end
        
        function J = cost(obj, x_all, u_all)
           J = 0;
           for itr =  1 : obj.sim_N-1
               J = J + x_all(itr, :)*obj.Q*x_all(itr, :)' + u_all(itr, :)*obj.R*u_all(itr, :)';
           end
        end
        
        function R = reward(obj, x, u)
            lidx = 2; % rectifier の周波数成分のみ小さくすれば良い
            R = -1/10*(x(lidx)*x(lidx) + u*obj.R*u');
        end
        
        function [local_x_all, env_x_all, rect_x_all, y_xhat_w_v_all, rl_u_all] = sim(obj, theta_mu, seed)
            if nargin < 3
                seed = rng();
            end
            rng(seed);
            if nargin < 2 || isempty(theta_mu)
                theta_mu = obj.policy.get_params();
            end
            d_L = randn(obj.sim_N, 2);
            local_x_all = zeros(obj.sim_N, obj.model.local_nx);
            env_x_all = zeros(obj.sim_N, obj.model.env_nx);
            rect_x_all = zeros(obj.sim_N, obj.model.rect_nx);
            y_xhat_w_v_all = zeros(obj.sim_N, obj.model.ny+obj.model.ny+obj.model.nw+obj.model.nv);
            rl_u_all = zeros(obj.sim_N, obj.model.nu);
            mpc_u_all = zeros(obj.sim_N, obj.model.nu);
            % calculate MBC gain
            if obj.model.apx_nx == 0
%                 K = lqr(obj.model.A, obj.model.B, obj.Q, obj.R);
                obj.model.set_controlled_system(obj.Q, obj.R);
            else
%                 K = lqr(obj.model.A, obj.model.B, blkdiag(obj.Q, eye(obj.model.apx_nx)*obj.Q1), obj.R);
                obj.model.set_controlled_system(blkdiag(obj.Q, eye(obj.model.apx_nx)*obj.Q1), obj.R);
            end
%             K1 = K(1:obj.model.local_nx);
%             K2 = K(obj.model.local_nx+1 : end);
%             % set episode initial
%             local_x_all(1, :) = ini';% When network model, local system state
            for k = 1 : obj.sim_N-1
               % MBC input
               mpc_u_all(k, :) = 0; % LQR Controller already implemented
%                    mpc_u_all(k, :) = ([K1,-K2]*rect_x_all(k, :)')' -(K1*ywv_all(k, 1:obj.model.ny)')';
               % RL input
               rl_u_all(k, :) = obj.policy.determistic_policy([local_x_all(k, :), rect_x_all(k, :)], theta_mu);
               %  observe S_(k+1)
               [local_ne_x, env_ne_x, ne_ywv, rect_ne_x] = ...
                    obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)',rect_x_all(k, :)', y_xhat_w_v_all(k, :)', rl_u_all(k, :) + mpc_u_all(k, :), d_L(k, :)');

                local_x_all(k+1, :) = local_ne_x';
               env_x_all(k+1, :) = env_ne_x';
               y_xhat_w_v_all(k+1, :) = ne_ywv';
               rect_x_all(k+1, :) = rect_ne_x';
           end
        end
%         
        function [local_x_all, env_x_all, rect_x_all, y_xhat_w_v_all] = sim_lqrcontroller(obj, seed)
            if nargin < 2
                seed = rng();
            end
            rng(seed);
            d_L = randn(obj.sim_N, 2);
            local_x_all = zeros(obj.sim_N, obj.model.local_nx);
            env_x_all = zeros(obj.sim_N, obj.model.env_nx);
            rect_x_all = zeros(obj.sim_N, obj.model.rect_nx);
            y_xhat_w_v_all = zeros(obj.sim_N, obj.model.ny+obj.model.ny+obj.model.nw+obj.model.nv);
            % calculate MBC gain
            if obj.model.apx_nx == 0
%                 K = lqr(obj.model.A, obj.model.B, obj.Q, obj.R);
                obj.model.set_controlled_system(obj.Q, obj.R);
            else
%                 K = lqr(obj.model.A, obj.model.B, blkdiag(obj.Q, eye(obj.model.apx_nx)*obj.Q1), obj.R);
                obj.model.set_controlled_system(blkdiag(obj.Q, eye(obj.model.apx_nx)*obj.Q1), obj.R);
            end
%             K1 = K(1:obj.model.local_nx);
%             K2 = K(obj.model.local_nx+1 : end);
%             % set episode initial
%             local_x_all(1, :) = ini';% When network model, local system state
            for k = 1 : obj.sim_N-1
               %  observe S_(k+1)
               [local_ne_x, env_ne_x, ne_ywv, rect_ne_x] = ...
                    obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)',rect_x_all(k, :)', y_xhat_w_v_all(k, :)', 0, d_L(k, :)');

                local_x_all(k+1, :) = local_ne_x';
               env_x_all(k+1, :) = env_ne_x';
               y_xhat_w_v_all(k+1, :) = ne_ywv';
               rect_x_all(k+1, :) = rect_ne_x';
           end 
        end
        
        % netwerk だと生のvwの計測データがないと初期点すら計算できない．
%         
%         function w = set_w(obj, K)
%             if nargin < 3
%                 w0 = zeros(obj.value.apx_function.N, 1);
%             end
%             [X,Y] = meshgrid(-1:.1:1, -1:.1:1);
%             mesh_size = size(X, 1);
%             Z = zeros(mesh_size, mesh_size);
%             for itr1 = 1 : mesh_size
%                for itr2 = 1 : mesh_size
%                     Z(itr1, itr2) = obj.cumulative_reward(K, [X(itr1,itr2);Y(itr1,itr2)]);
%                end
%             end
%             state = [vec(X), vec(Y)];
%             target = vec(Z);
% 
%             options = optimoptions('lsqnonlin','Display','iter','SpecifyObjectiveGradient',true);
%             w = lsqnonlin(@(w)obj.apx_cost_function(state, target, w),w0,[],[],options);
%         end
%         
%         function [r ,apx_all, u_all] = cumulative_reward(obj, K, ini)
%             apx_all = zeros(obj.sim_N, obj.model.apx_nx);
%             y_all = zeros(obj.sim_N, obj.model.ny);
%             u_all = zeros(obj.sim_N, obj.model.nu);
%             r = 0;
%             % set initial
%             apx_all(1, :) = ini';
%             for k = 1 : obj.sim_N-1
%                 u_all(k, :) = (-K*apx_all(k, :)')';
%                 [ne_x, y] = obj.model.apx_dynamics(apx_all(k, :)', u_all(k, :)');
%                 apx_all(k+1, :) = ne_x';
%                 y_all(k, :) = y';
%                 r = r + obj.gamma^(k-1)*obj.reward(apx_all(k+1, :), u_all(k, :));
%             end
%         end
%         
%         function [error, grad] = apx_cost_function(obj, x, target, w)
%             apx = zeros(size(target));
%             for itr = 1 : length(target)
%                 apx(itr) = obj.value.est_value(x(itr, :), w);
%             end
%             error = apx - target;
%             if nargout > 1
%                 grad = zeros(length(target), obj.value.apx_function.N);
%                 for itr = 1 : length(target)
%                     grad(itr, :) = obj.value.value_grad(x(itr, :));
%                 end
%             end
%         end
    end
end

