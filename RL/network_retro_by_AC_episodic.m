classdef network_retro_by_AC_episodic < RL_train
    % レトロフィット制御問題に対するRL_actor_criticのMore General class
    
    properties
        opt_policy
        opt_value
        belief_sys
        Q1 = 1
        Q_c = diag([1,1])
        Q_r = diag([1,1])
        R = 1
        gamma = 0.99
        max_episode = 3.5e3
        snapshot = 100;
    end
    
    methods
        function obj = network_retro_by_AC_episodic(model, opt_policy, opt_value, belief_sys)
            obj.model = model;
            obj.opt_policy = opt_policy;
            obj.opt_value  = opt_value;
            obj.belief_sys = belief_sys;
            obj.Q1 = 1;
            obj.Q_c = diag([1,1]);
            obj.Q_r = diag([1,1]);
            obj.R = 1;
            obj.gamma = 0.99;
            obj.max_episode = 3.5e3;
            obj.snapshot = 100;
        end
        
        function [local_x_all, mpc_u_all, rl_u_all, policy_snapshot, value_snapshot, reward_history, varargout] = train(obj, ini, Te, seed, varargin)
            sim_N = Te/obj.model.Ts + 1;
            t = (0:obj.model.Ts:Te)';
            if nargin < 3 || isempty(seed)
                seed = rng();
            end
            rng(seed)
            reward_history = zeros(obj.max_episode, 1);
            % record point
            record_point = obj.snapshot:obj.snapshot:obj.max_episode;
            value_snapshot = zeros(obj.opt_value.approximate_function_class.apx_function.N, length(record_point));
            policy_snapshot = zeros(obj.opt_policy.approximate_function_class.apx_function.N, length(record_point));
            % calculate MBC gain
            K = zeros(size(obj.model.B'));
            tmp = strcmp(varargin, 'parallel');
            if sum(tmp)
                if  varargin{find(tmp)+1}
                    K = lqr(obj.model.A, obj.model.B, blkdiag(obj.Q_c, eye(obj.model.apx_nx)*obj.Q1), obj.R);
                end
            end
            K1 = K(1:obj.model.local_nx);
            K2 = K(obj.model.local_nx+1 : end);
            clearvars tmp
            % local noise
            d_L = zeros(sim_N, 2);
            % generate initial
            local_ini = repmat(ini, obj.max_episode, 1);
%             local_ini = rand(obj.max_episode, obj.model.local_nx) - 0.5;
            env_ini   = zeros(obj.max_episode, obj.model.env_nx);
            rect_ini  = zeros(obj.max_episode, obj.model.rect_nx);
            % start episode learning
            record_idx = 1;
            belief_state = zeros(2, size(obj.belief_sys, 1));% belief state % 1 line: current state , 2 line: previous state
            for episode = 1 : obj.max_episode
                % memory reset
                local_x_all = nan(sim_N, obj.model.local_nx);
                env_x_all = nan(sim_N, obj.model.env_nx);
                rect_x_all = nan(sim_N, obj.model.rect_nx);
                y_w_v_all = nan(sim_N, obj.model.ny+obj.model.nw+obj.model.nv);
                rl_u_all = nan(sim_N, obj.model.nu);
                mpc_u_all = nan(sim_N, obj.model.nu);
                % episode initialize
                obj.opt_policy.initialize();
                obj.opt_value.initialize();
                reward = 0;
                % belief initialize
                belief_state = zeros(size(belief_state));
                % set episode initial
                local_x_all(1, :) = local_ini(episode, :);
                env_x_all(1, :) = env_ini(episode, :);
                rect_x_all(1, :) = rect_ini(episode, :);
                % simlate
                for k = 1 : sim_N-1
                    % MBC input
                    % LQR Controller 
                    ywv = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)');
                    y_w_v_all(k, :) = ywv';
                    mpc_u_all(k, :) = ([K1,-K2]*rect_x_all(k, :)')' -(K1*y_w_v_all(k, 1:obj.model.ny)')';
                    % RL input
                    rl_u_all(k, :) = obj.opt_policy.approximate_function_class.stocastic_policy(belief_state(1, :), []);
                    %  observe State(k+1)
                    [~, local_ne_x, env_ne_x] = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)', rl_u_all(k, :) + mpc_u_all(k, :), d_L(k, :)');
                    local_x_all(k+1, :) = local_ne_x';
                    env_x_all(k+1, :) = env_ne_x';
                    [rect_ne_x, rect_yw] = obj.model.rect_dynamics(local_x_all(k, :)', env_x_all(k, :)', rect_x_all(k, :)');
                    rect_x_all(k+1, :) = rect_ne_x';
                    % belief upadate
                    belief_state(2, :) = belief_state(1, :);
                    belief_state(1, :) = (obj.belief_sys*[belief_state(1, :)'; rect_yw])'; % momory store
                    % Get Reward r
                    r = obj.reward([local_x_all(k+1, :), rect_x_all(k+1, :)], rl_u_all(k, :)+mpc_u_all(k, :));
                    reward =  reward + obj.gamma^(k-1)*r;
                    % TD Erorr
                    V_k1 = obj.opt_value.approximate_function_class.est_value(belief_state(1, :)); 
                    V_k0 = obj.opt_value.approximate_function_class.est_value(belief_state(2, :));
                    delta = r + obj.gamma*V_k1 - V_k0;
                    % parameter update
                    obj.opt_policy.opt(delta, belief_state(2, :), rl_u_all(k, :), obj.gamma, obj.model, obj.belief_sys);
                    obj.opt_value.opt(delta, belief_state(2, :), rl_u_all(k, :), obj.gamma);
                end
                % record history
                if ~mod(episode, obj.snapshot)
                    value_snapshot(:, record_idx) = obj.opt_value.approximate_function_class.get_params(); 
                    policy_snapshot(:, record_idx) = obj.opt_policy.approximate_function_class.get_params();
                    record_idx =  record_idx + 1;
                end
                reward_history(episode) = reward;
                % callback
                subplot(3,1,1)
                plot(t, local_x_all)
                title(['\fontsize{16}','Episode-',num2str(episode)])
                ylabel('z')
                grid on
                xlim([0, t(end)])
                lim = max(max(abs(local_x_all),[],'omitnan'));
                if isempty(lim) || lim == 0
                    lim = 1;
                end
                ylim([-lim, lim]);
                subplot(3,1,2)
                plot(t, rect_x_all)
                grid on
                ylabel('\hat{z}')
                xlim([0, t(end)])
                lim = max(max(abs(rect_x_all),[],'omitnan'));
                if isempty(lim) || lim == 0 || isnan(lim)
                    lim = 1;
                end
                ylim([-lim, lim]);
                subplot(3,1,3)
                plot(nonzeros(reward_history),'-b')
                ylabel('Culumative Reward')
                drawnow
                disp(strcat('Episode-',num2str(episode),' : value  constraint update times : ', num2str(obj.opt_value.counter) ,'/',num2str(sim_N-1)))
                disp(strcat('Episode-',num2str(episode),' : policy constraint update times : ', num2str(obj.opt_policy.counter) ,'/',num2str(sim_N-1)))
                if nargout > 7
                   varargout{1}(episode) = getframe(gcf);
               end
            end
        end
        
        function J = cost(obj, x_all, u_all)
           J = 0;
           for itr =  1 : sim_N-1
               J = J + x_all(itr, :)*obj.Q_c*x_all(itr, :)' + u_all(itr, :)*obj.R*u_all(itr, :)';
           end
        end
        
        function R = reward(obj, x, u)
            lidx = 2; % rectifier の周波数成分のみ小さくすれば良い
%             R = -1/10*(x(lidx)*obj.Q_r(2,2)*x(lidx) + u*obj.R*u');
%             R = -1/10*(x*blkdiag(obj.Q_r, obj.Q_r, eye(obj.model.apx_nx)*obj.Q1)*x' + u*obj.R*u');
            R = -1/10*(x(1:2)*obj.Q_r*x(1:2)' + u*obj.R*u');
        end
        
        function [local_x_all, env_x_all, rect_x_all, y_w_v_all, rl_u_all, reward] = sim(obj, ini, Te, theta_mu, noise_power, seed, varargin)
            if nargin < 6
                seed = rng();
            end
            rng(seed);
            if nargin < 3 || isempty(theta_mu)
                theta_mu = obj.opt_policy.approximate_function_class.get_params();
            end
            if nargin < 2 || isempty(ini)
                ini = randn(obj.model.local_nx, 1);
            end
            if nargin < 5 || isempty(noise_power)
                noise_power = 1;
            end
            sim_N = Te/obj.model.Ts + 1;
            t = (0:obj.model.Ts:Te)';
            d_L = noise_power*randn(sim_N, 2);
            local_x_all = zeros(sim_N, obj.model.local_nx);
            env_x_all = zeros(sim_N, obj.model.env_nx);
            rect_x_all = zeros(sim_N, obj.model.rect_nx);
            y_w_v_all = zeros(sim_N, obj.model.ny+obj.model.nw+obj.model.nv);
            rl_u_all = zeros(sim_N, obj.model.nu);
            mpc_u_all = zeros(sim_N, obj.model.nu);
            obj.opt_policy.approximate_function_class.initialize_memory();
            % set initial
            local_x_all(1, :) = ini';
            % calculate MBC gain
            K = lqr(obj.model.A, obj.model.B, blkdiag(obj.Q_c, eye(obj.model.apx_nx)*obj.Q1), obj.R);
            tmp = strcmp(varargin, 'parallel');
            if sum(tmp)
                if  strcmp(varargin{find(tmp)+1},'off')
                    K = zeros(size(K));
                end
            end
            K1 = K(1:obj.model.local_nx);
            K2 = K(obj.model.local_nx+1 : end);
            % observe belief state
            belief_state = zeros(1, size(obj.belief_sys, 1));% belief state % 1 line: current state , 2 line: previous state
            % initialize reward
            reward = 0;
            for k = 1 : sim_N-1
                % LQR Controller 
                ywv = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)');
                y_w_v_all(k, :) = ywv';
                mpc_u_all(k, :) = ([K1,-K2]*rect_x_all(k, :)')' -(K1*y_w_v_all(k, 1:obj.model.ny)')';
                % RL input
                rl_u_all(k, :) = obj.opt_policy.approximate_function_class.stocastic_policy(belief_state(1, :), theta_mu);
                %  observe S_(k+1)
                [~, local_ne_x, env_ne_x] = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)', rl_u_all(k, :) + mpc_u_all(k, :), d_L(k, :)');
                local_x_all(k+1, :) = local_ne_x';
                env_x_all(k+1, :) = env_ne_x';
                [rect_ne_x, rect_yw] = obj.model.rect_dynamics(local_x_all(k, :)', env_x_all(k, :)', rect_x_all(k, :)');
                rect_x_all(k+1, :) = rect_ne_x';
                % belief_update
                belief_state(1, obj.model.ny+1:end) = belief_state(1, 1:end-obj.model.ny);
                belief_state(1, 1:obj.model.ny) = rect_yw'; % momory store
                r = obj.reward([local_x_all(k+1, :), rect_x_all(k+1, :)], rl_u_all(k, :)+mpc_u_all(k, :));
                reward =  reward + obj.gamma^(k-1)*r;
           end
        end
%         
        function [local_x_all, env_x_all, rect_x_all, y_w_v_all, reward] = sim_lqrcontroller(obj, ini, Te, noise_power, seed)
            if nargin < 4
                seed = rng();
            end
            if nargin < 3 ||isempty(noise_power)
                noise_power = 1;
            end
            sim_N = Te/obj.model.Ts + 1;
            t = (0:obj.model.Ts:Te)';
            rng(seed);
            d_L = noise_power*randn(sim_N, 2);
            local_x_all = zeros(sim_N, obj.model.local_nx);
            env_x_all = zeros(sim_N, obj.model.env_nx);
            rect_x_all = zeros(sim_N, obj.model.rect_nx);
            y_w_v_all = zeros(sim_N, obj.model.ny+obj.model.nw+obj.model.nv);
            mpc_u_all = zeros(sim_N, obj.model.nu);
            % set initial
            local_x_all(1, :) = ini';
            % calculate MBC gain
            K = lqr(obj.model.A, obj.model.B, blkdiag(obj.Q_c, eye(obj.model.apx_nx)*obj.Q1), obj.R);
            K1 = K(1:obj.model.local_nx);
            K2 = K(obj.model.local_nx+1 : end);
            % initialize reward
            reward = 0;
            for k = 1 : sim_N-1
                 % LQR Controller 
                ywv = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)');
                y_w_v_all(k, :) = ywv';
                mpc_u_all(k, :) = ([K1,-K2]*rect_x_all(k, :)')' -(K1*y_w_v_all(k, 1:obj.model.ny)')';
                %  observe S_(k+1)
                [~, local_ne_x, env_ne_x] = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)', mpc_u_all(k, :), d_L(k, :)');
                local_x_all(k+1, :) = local_ne_x';
                env_x_all(k+1, :) = env_ne_x';
                [rect_ne_x, rect_yw] = obj.model.rect_dynamics(local_x_all(k, :)', env_x_all(k, :)', rect_x_all(k, :)');
                rect_x_all(k+1, :) = rect_ne_x';
                r = obj.reward([local_x_all(k+1, :), rect_x_all(k+1, :)], mpc_u_all(k, :));
                reward =  reward + obj.gamma^(k-1)*r;
            end 
        end
        
        function [local_x_all, env_x_all, y_w_v_all, reward] = sim_original(obj, ini, Te, noise_power, seed)
            if nargin < 4
                seed = rng();
            end
            if nargin < 3 || isempty(noise_power)
                noise_power = 1;
            end
            sim_N = Te/obj.model.Ts + 1;
            t = (0:obj.model.Ts:Te)';
            rng(seed);
            d_L = noise_power*randn(sim_N, 2);
            local_x_all = zeros(sim_N, obj.model.local_nx);
            env_x_all = zeros(sim_N, obj.model.env_nx);
            y_w_v_all = zeros(sim_N, obj.model.ny+obj.model.nw+obj.model.nv);
            % set initial
            local_x_all(1, :) = ini';
            % initialize reward
            reward = 0;
            for k = 1 : sim_N-1
                ywv = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)');
                y_w_v_all(k, :) = ywv';
                %  observe S_(k+1)
                [~, local_ne_x, env_ne_x] = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)', 0, d_L(k, :)');
                local_x_all(k+1, :) = local_ne_x';
                env_x_all(k+1, :) = env_ne_x';
                r = obj.reward(local_x_all(k+1, :), 0);
                reward =  reward + obj.gamma^(k-1)*r;
            end 
        end
    end
end

