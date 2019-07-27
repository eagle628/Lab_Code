classdef netwrok_retro_by_adaptive_and_POMDP_episodic < RL_train
    % レトロフィ�?ト制御問題に対するRL_adaptive_LQR_RLSのclass
    
    properties(Constant)        
        gamma = 1
        max_episode = 2000%3.5e3
        snapshot = 100;
    end
    
    properties
        sim_N
        t
        belief_N
        Q
        R
    end
    
    methods
        function obj = netwrok_retro_by_adaptive_and_POMDP_episodic(model, belief_N, Te)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            obj.belief_N = belief_N;
            obj.Q = eye(obj.belief_N*obj.model.ny);
            obj.R = eye(obj.model.nu);
        end
        
        function [local_x_all, rl_u_all, theta_snapshot, reward_history] = train(obj, ini, seed)
            if nargin < 3
                seed = rng();
            end
            rng(seed)
            reward_history = zeros(obj.max_episode, 1);
            local_x_all = zeros(obj.sim_N, obj.model.local_nx);
            env_x_all = zeros(obj.sim_N, obj.model.env_nx);
            rect_x_all = zeros(obj.sim_N, obj.model.rect_nx);
            y_w_v_all = zeros(obj.sim_N, obj.model.ny+obj.model.nw+obj.model.nv);
            rl_u_all = zeros(obj.sim_N, obj.model.nu);
            % set episode initial
            local_x_all(1, :) = ini';% When network model, local system state
            % local noise
            d_L = zeros(obj.sim_N, 2);
            % parma initialize
            n_th = (obj.model.ny*obj.belief_N+obj.model.nu)*(obj.model.ny*obj.belief_N+obj.model.nu+1)/2;
            theta = zeros(n_th, 1); 
            % record point
            record_point = obj.snapshot:obj.snapshot:obj.max_episode;
            theta_snapshot = zeros(size(theta, 1), length(record_point));
            % start episode learning
            record_idx = 1;
            belief_state = zeros(2, obj.model.ny*obj.belief_N);% belief state % 1 line: current state , 2 line: previous state
            P = eye(n_th)*10;
            K = zeros(obj.model.nu, obj.model.ny*obj.belief_N);
            sigma_pi = 1;
            for episode = 1 : obj.max_episode
                reward = 0;
                belief_state = zeros(size(belief_state));
                for k = 1 : obj.sim_N-1
                    belief_state(2, :) = belief_state(1, :);
                    belief_state(1, obj.model.ny+1:end) = belief_state(1, 1:end-obj.model.ny);
                    % current observe
                    ywv = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)');
                    y_w_v_all(k, :) = ywv';
%                     mpc_u_all(k, :) = ([K1,-K2]*rect_x_all(k, :)')' -(K1*y_w_v_all(k, 1:obj.model.ny)')';
                    % RL input
                    rl_u_all(k, :) = K*belief_state(1, :)' + sigma_pi*randn(obj.model.nu, 1);
                    %  observe S_(k+1)
                    [~, local_ne_x, env_ne_x] = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)', rl_u_all(k, :), d_L(k, :)');
                    local_x_all(k+1, :) = local_ne_x';
                    env_x_all(k+1, :) = env_ne_x';
                    [rect_ne_x, rect_yw] = obj.model.rect_dynamics(local_x_all(k, :)', env_x_all(k, :)', rect_x_all(k, :)');
                    rect_x_all(k+1, :) = rect_ne_x';
                    belief_state(1, 1:obj.model.ny) = rect_yw(1:obj.model.ny)';
                    % Get Reward r
                    r = obj.reward(belief_state(1, :)', rl_u_all(k, :));
                    reward =  reward + obj.gamma^(k-1)*r;
                    if k ~= 1 && episode ~= 1
                        % TD Erorr for Recruisive least square
                        phi = H_to_theta([belief_state(2, :)';rl_u_all(k-1, :)';]*[belief_state(2, :)';rl_u_all(k-1, :)';]')...
                                - obj.gamma*H_to_theta([belief_state(1, :)';K*belief_state(1, :)';]*[belief_state(1, :)';K*belief_state(1, :)';]');
                        epsilon = (-r) - phi' * theta;
                        denom = phi'*P*phi;
                        theta = theta + P*phi*epsilon/denom;
                        P = P - (P*phi*(phi')*P)/denom;
                    end
                end
                if episode ~= 1
                    H = theta_to_H(theta, obj.model.ny*obj.belief_N+obj.model.nu);
                    K = - inv(H(obj.model.ny*obj.belief_N+1:obj.model.ny*obj.belief_N+obj.model.nu,obj.model.ny*obj.belief_N+1:obj.model.ny*obj.belief_N+obj.model.nu))...
                            *H(obj.model.ny*obj.belief_N+1:obj.model.ny*obj.belief_N+obj.model.nu,1:obj.model.ny*obj.belief_N);
                end
                % record history
                if ~mode(episode, obj.snapshot)
                    theta_snapshot(:, record_idx) = theta;
                    record_idx =  record_idx + 1;
                end
                reward_history(episode) = reward;
                % callback
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
%                 ylim([-5,0])
                drawnow
%                 cost_history(episode) = obj.cost([local_x_all(k, :), rect_x_all(k, :)], rl_u_all+mpc_u_all);
%                 figure(1)
%                 callback_RL(episode, obj.t, [local_x_all(k, :), rect_x_all(k, :)], cost_history, reward_history)
            end
        end
        
        function J = cost(obj, x_all, u_all)
           J = 0;
           for itr =  1 : obj.sim_N-1
               J = J + x_all(itr, :)*obj.Q_c*x_all(itr, :)' + u_all(itr, :)*obj.R*u_all(itr, :)';
           end
        end
        
        function R = reward(obj, x, u)
            R = -(x'*obj.Q*x + u'*obj.R*u);
        end
        
        function [local_x_all, env_x_all, rect_x_all, y_xhat_w_v_all, rl_u_all] = sim(obj, theta_mu, ini, noise_power, seed)
            if nargin < 5
                seed = rng();
            end
            rng(seed);
            if nargin < 2 || isempty(theta_mu)
                theta_mu = obj.policy.get_params();
            end
            d_L = noise_power*randn(obj.sim_N, 2);
            local_x_all = zeros(obj.sim_N, obj.model.local_nx);
            env_x_all = zeros(obj.sim_N, obj.model.env_nx);
            rect_x_all = zeros(obj.sim_N, obj.model.rect_nx);
            y_xhat_w_v_all = zeros(obj.sim_N, obj.model.ny+obj.model.ny+obj.model.nw+obj.model.nv);
            rl_u_all = zeros(obj.sim_N, obj.model.nu);
            mpc_u_all = zeros(obj.sim_N, obj.model.nu);
            % set initial
            local_x_all(1, :) = ini';
            % calculate MBC gain
%                 K = lqr(obj.model.A, obj.model.B, blkdiag(obj.Q_c, eye(obj.model.apx_nx)*obj.Q1), obj.R);
            obj.model.set_controlled_system(blkdiag(obj.Q_c, eye(obj.model.apx_nx)*obj.Q1), obj.R);
%             K1 = K(1:obj.model.local_nx);
%             K2 = K(obj.model.local_nx+1 : end);
%             % set episode initial
%             local_x_all(1, :) = ini';% When network model, local system state
            for k = 1 : obj.sim_N-1
               % MBC input
               mpc_u_all(k, :) = 0; % LQR Controller already implemented
%                    mpc_u_all(k, :) = ([K1,-K2]*rect_x_all(k, :)')' -(K1*ywv_all(k, 1:obj.model.ny)')';
               % RL input
%                rl_u_all(k, :) = obj.policy.determistic_policy([local_x_all(k, :), rect_x_all(k, :)], theta_mu);
               rl_u_all(k, :) = obj.policy.stocastic_policy([local_x_all(k, :), rect_x_all(k, :)], theta_mu);
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
        function [local_x_all, env_x_all, rect_x_all, y_xhat_w_v_all, reward] = sim_lqrcontroller(obj, ini, noise_power, seed)
            if nargin < 4
                seed = rng();
            end
            rng(seed);
            d_L = noise_power*randn(obj.sim_N, 2);
            local_x_all = zeros(obj.sim_N, obj.model.local_nx);
            env_x_all = zeros(obj.sim_N, obj.model.env_nx);
            rect_x_all = zeros(obj.sim_N, obj.model.rect_nx);
            y_xhat_w_v_all = zeros(obj.sim_N, obj.model.ny+obj.model.ny+obj.model.nw+obj.model.nv);
            % set initial
            local_x_all(1, :) = ini';
            % calculate MBC gain
            K = lqr(obj.model.A, obj.model.B, blkdiag(obj.Q_c, eye(obj.model.apx_nx)*obj.Q1), obj.R);
            obj.model.set_controlled_system(blkdiag(obj.Q_c, eye(obj.model.apx_nx)*obj.Q1), obj.R);
            K1 = K(1:obj.model.local_nx);
            K2 = K(obj.model.local_nx+1 : end);
%             % set episode initial
%             local_x_all(1, :) = ini';% When network model, local system state
            reward = 0;
            for k = 1 : obj.sim_N-1
                %  observe S_(k+1)
                [local_ne_x, env_ne_x, ne_ywv, rect_ne_x] = ...
                    obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)',rect_x_all(k, :)', y_xhat_w_v_all(k, :)', 0, d_L(k, :)');

                local_x_all(k+1, :) = local_ne_x';
                env_x_all(k+1, :) = env_ne_x';
                y_xhat_w_v_all(k+1, :) = ne_ywv';
                rect_x_all(k+1, :) = rect_ne_x';
                r = obj.reward([local_x_all(k+1, :), rect_x_all(k+1, :)], ([K1,-K2]*rect_x_all(k, :)')' -(K1*y_xhat_w_v_all(k, 1:obj.model.ny)')');
                reward =  reward + obj.gamma^(k-1)*r;
            end 
        end
        
        % netwerk �?と生�?�vwの計測�?ータがな�?と初期点すら計算できな�?�?
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

%% local
function theta = H_to_theta(H)

    n = size(H);
    theta = H(1,:)';
    for i = 2:n
        theta = [theta; H(i,i:n)'];
    end
end

function H = theta_to_H(theta,n)
    H = theta(1:n)';
    k = n;
    for i = 2:n
        H = [ H ; [ zeros(1,i-1) theta(k+1:k+n-i+1)' ] ];
        k = k+n-i+1;
    end
    H = (H+H')/2;
end
