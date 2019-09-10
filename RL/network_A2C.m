classdef network_A2C < RL_train
    % 一般的な制御問題に対するRL_actor_criticのclass
    
    properties(Constant)
        Q1 = 1
        Q_c = diag([1,1])
        Q_r = diag([1,1])
        R = 1
        alpha = 1e-6
        beta = 1e-6
        gamma = 0.99
        lambda = 0.3 
        max_episode = 4e4
        snapshot = 100
    end
    
    properties
        policy
        value
        sim_N
        t
        advantage_N
        belief_N
    end
    
    methods
        function obj = network_A2C(model, policy, value, Te, advantage_N, belief_N)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            obj.policy = policy;
            obj.value  =  value;
            if nargin < 5
                advantage_N = 5;
            end
            obj.advantage_N = advantage_N;
            obj.belief_N = belief_N;
        end
        
        function [local_x_all, mpc_u_all, rl_u_all, theta_mu_snapshot, w_snapshot, reward_history, varargout] = train(obj, ini, seed, varargin)
            if nargin < 3 || isempty(seed)
                seed = rng();
            end
            rng(seed)
            rng('shuffle')
            cost_history = zeros(obj.max_episode, 1);
            reward_history = zeros(obj.max_episode, 1);
            % record point
            record_point = obj.snapshot:obj.snapshot:obj.max_episode;
            w_snapshot = zeros(obj.value.apx_function.N, length(record_point));
            theta_mu_snapshot = zeros(obj.policy.apx_function.N, length(record_point));
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
            % local noise
            d_L = zeros(obj.sim_N, 2);
            % params initialize
            w = obj.value.get_params();
            theta_mu = obj.policy.get_params();
            % start episode learning
            record_idx = 1;
            belief_state = zeros(obj.advantage_N+1, obj.belief_N*obj.model.ny);% belief state % 1 line: current state , 2 line: previous state
            for episode = 1 : obj.max_episode
                reward = 0;
                % memory reset
                local_x_all = nan(obj.sim_N, obj.model.local_nx);
                env_x_all = nan(obj.sim_N, obj.model.env_nx);
                rect_x_all = nan(obj.sim_N, obj.model.rect_nx);
                y_w_v_all = nan(obj.sim_N, obj.model.ny+obj.model.nw+obj.model.nv);
                rl_u_all = nan(obj.sim_N, obj.model.nu);
                mpc_u_all = nan(obj.sim_N, obj.model.nu);
                constraint_update = 0;
                obj.policy.initialize_memory();
                % belief initialize
                belief_state = zeros(size(belief_state));
                % set episode initial
                local_x_all(1, :) = ini';% When network model, local system state
%                 local_x_all(1, :) = rand(1,2)-0.5;
                env_x_all(1, :) = zeros(1, obj.model.env_nx);
                rect_x_all(1, :) = zeros(1, obj.model.rect_nx);
                for k = 1 : obj.sim_N-1
                    % MBC input
                    % LQR Controller 
                    ywv = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)');
                    y_w_v_all(k, :) = ywv';
                    mpc_u_all(k, :) = ([K1,-K2]*rect_x_all(k, :)')' -(K1*y_w_v_all(k, 1:obj.model.ny)')';
                    % RL input
                    tmp = strcmp(varargin, 'Input-Clipping');
                    if sum(tmp)
                        rl_u_all(k, :) = obj.policy.stocastic_policy(belief_state(1, :), [], 'Input-Clipping', varargin{find(tmp)+1});
                    else
                        rl_u_all(k, :) = obj.policy.stocastic_policy(belief_state(1, :), []);
                    end
                    %  observe State(k+1)
                    [~, local_ne_x, env_ne_x] = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)', rl_u_all(k, :) + mpc_u_all(k, :), d_L(k, :)');
                    local_x_all(k+1, :) = local_ne_x';
                    env_x_all(k+1, :) = env_ne_x';
                    [rect_ne_x, rect_yw] = obj.model.rect_dynamics(local_x_all(k, :)', env_x_all(k, :)', rect_x_all(k, :)');
                    rect_x_all(k+1, :) = rect_ne_x';
                    % belief upadate
                    belief_state(2:end, :) = belief_state(1:end-1, :);
                    belief_state(1, obj.model.ny+1:end) = belief_state(1, 1:end-obj.model.ny);
                    belief_state(1, 1:obj.model.ny) = rect_yw'; % momory store
                    if ~mod(k, obj.advantage_N)
                        % reset accumulate gradients
                        d_theta_mu = zeros(size(theta_mu));
                        d_w        = zeros(size(w));
                        backward_reward = obj.value.est_value(belief_state(1, :));
                        belief_idx = 1;
                        for iter1 = k : -1 : k - obj.advantage_N+1
                            % Get Reward r
                            r = obj.reward(belief_state(belief_idx, :), rl_u_all(iter1, :));
                            backward_reward = r + obj.gamma*backward_reward;
                            advantage = backward_reward - obj.value.est_value(belief_state(belief_idx+1, :)); 
                            % accumulates gradients
                            d_theta_mu = d_theta_mu + obj.policy.policy_grad_mu(rl_u_all(belief_idx+1, :), belief_state(belief_idx+1, :))*advantage;
                            d_w        = d_w        + 2*obj.value.value_grad(belief_state(belief_idx+1, :))*advantage;
% % %                             % Get Reward r
% % %                             r = obj.reward(x_all(iter1+1, :), rl_u_all(iter1, :));
% % %                             backward_reward = r + obj.gamma*backward_reward;
% % %                             backward_reward_lambda = r + obj.gamma*(obj.lambda*backward_reward+(1-obj.lambda)*obj.value.est_value(x_all(iter1, :)));
% % %                             % accumulates gradients
% % %                             d_theta_mu = d_theta_mu + obj.policy.policy_grad_mu(rl_u_all(iter1, :), x_all(iter1, :))*(backward_reward - obj.value.est_value(x_all(iter1, :)));
% % %                             d_w        = d_w        - 2*obj.value.value_grad(x_all(iter1, :))*(backward_reward_lambda - obj.value.est_value(x_all(iter1, :)));
% % %                             % update params
% % %                             obj.value.set_params(w+obj.alpha*d_w);
% % %                             obj.policy.set_params(theta_mu+obj.beta*d_theta_mu);
% % %                             w = obj.value.get_params();
% % %                             theta_mu = obj.policy.get_params();
                            belief_idx = belief_idx + 1;
                        end
                        % update params
                        tmp = strcmp(varargin, 'Invalid-Constraint');
                        if obj.policy.policy_constraint(theta_mu+obj.beta*d_theta_mu, obj.model, obj.belief_N, 'Invalid-Constraint', varargin{find(tmp)+1})
                            obj.value.set_params(w+obj.alpha*d_w);
                            obj.policy.set_params(theta_mu+obj.beta*d_theta_mu);
                            w = obj.value.get_params();
                            theta_mu = obj.policy.get_params();
                            constraint_update = constraint_update + 1;
                        end
                    end
                end
%                close 3
                % record history
                if ~mod(episode, obj.snapshot) 
                    w_snapshot(:, record_idx) = w; 
                    theta_mu_snapshot(:, record_idx) = theta_mu;
                    record_idx =  record_idx + 1;
                end
                reward_history(episode) = reward;
                figure(1)
                % callback
                subplot(3,1,1)
                plot(obj.t, local_x_all)
                title(['\fontsize{16}','Episode-',num2str(episode)])
                ylabel('z')
                grid on
                xlim([0, obj.t(end)])
                lim = max(max(abs(local_x_all),[],'omitnan'));
                if isempty(lim) || lim == 0
                    lim = 1;
                end
                ylim([-lim, lim]);
                subplot(3,1,2)
                plot(obj.t, rect_x_all)
                grid on
                ylabel('\hat{z}')
                xlim([0, obj.t(end)])
                lim = max(max(abs(rect_x_all),[],'omitnan'));
                if isempty(lim) || lim == 0 || isnan(lim)
                    lim = 1;
                end
                ylim([-lim, lim]);
                subplot(3,1,3)
                plot(nonzeros(reward_history),'-b')
                ylabel('Culumative Reward')
                drawnow
                disp(strcat('Episode-',num2str(episode),' : update times : ', num2str(constraint_update) ,'/',num2str(obj.sim_N-1)))
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
%             R = -1/10*(x(lidx)*obj.Q_r(2,2)*x(lidx) + u*obj.R*u');
%             R = -1/10*(x*blkdiag(obj.Q_r, obj.Q_r, eye(obj.model.apx_nx)*obj.Q1)*x' + u*obj.R*u');
            R = -1/10*(x(1:2)*obj.Q_r*x(1:2)' + u*obj.R*u');
        end
        
        function [local_x_all, env_x_all, rect_x_all, y_w_v_all, rl_u_all, reward] = sim(obj, ini, theta_mu, noise_power, seed, varargin)
            if nargin < 6
                seed = rng();
            end
            rng(seed);
            if nargin < 3 || isempty(theta_mu)
                theta_mu = obj.policy.get_params();
            end
            if nargin < 2 || isempty(ini)
                ini = zeros(obj.model.local_nx, 1);
            end
            if nargin < 5 || isempty(noise_power)
                noise_power = 1;
            end
            d_L = noise_power*randn(obj.sim_N, 2);
            local_x_all = zeros(obj.sim_N, obj.model.local_nx);
            env_x_all = zeros(obj.sim_N, obj.model.env_nx);
            rect_x_all = zeros(obj.sim_N, obj.model.rect_nx);
            y_w_v_all = zeros(obj.sim_N, obj.model.ny+obj.model.nw+obj.model.nv);
            rl_u_all = zeros(obj.sim_N, obj.model.nu);
            mpc_u_all = zeros(obj.sim_N, obj.model.nu);
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
            belief_state = zeros(1, obj.belief_N*obj.model.ny);% belief state % 1 line: current state , 2 line: previous state
            % initialize reward
            reward = 0;
            for k = 1 : obj.sim_N-1
                % LQR Controller 
                ywv = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)');
                y_w_v_all(k, :) = ywv';
                mpc_u_all(k, :) = ([K1,-K2]*rect_x_all(k, :)')' -(K1*y_w_v_all(k, 1:obj.model.ny)')';
                % RL input
                rl_u_all(k, :) = obj.policy.stocastic_policy(belief_state(1, :), theta_mu);
                %  observe S_(k+1)
                [~, local_ne_x, env_ne_x] = obj.model.dynamics(local_x_all(k, :)', env_x_all(k, :)', rl_u_all(k, :) + mpc_u_all(k, :), d_L(k, :)');
                local_x_all(k+1, :) = local_ne_x';
                env_x_all(k+1, :) = env_ne_x';
                [rect_ne_x, rect_yw] = obj.model.rect_dynamics(local_x_all(k, :)', env_x_all(k, :)', rect_x_all(k, :)');
                rect_x_all(k+1, :) = rect_ne_x';
                belief_state(1, obj.model.ny+1:end) = belief_state(1, 1:end-obj.model.ny);
                belief_state(1, 1:obj.model.ny) = rect_yw'; % momory store
                r = obj.reward([local_x_all(k+1, :), rect_x_all(k+1, :)], rl_u_all(k, :)+mpc_u_all(k, :));
                reward =  reward + obj.gamma^(k-1)*r;
           end
        end
%         
        function [local_x_all, env_x_all, rect_x_all, y_w_v_all, reward] = sim_lqrcontroller(obj, ini, noise_power, seed)
            if nargin < 4
                seed = rng();
            end
            if nargin < 3 ||isempty(noise_power)
                noise_power = 1;
            end
            rng(seed);
            d_L = noise_power*randn(obj.sim_N, 2);
            local_x_all = zeros(obj.sim_N, obj.model.local_nx);
            env_x_all = zeros(obj.sim_N, obj.model.env_nx);
            rect_x_all = zeros(obj.sim_N, obj.model.rect_nx);
            y_w_v_all = zeros(obj.sim_N, obj.model.ny+obj.model.nw+obj.model.nv);
            mpc_u_all = zeros(obj.sim_N, obj.model.nu);
            % set initial
            local_x_all(1, :) = ini';
            % calculate MBC gain
            K = lqr(obj.model.A, obj.model.B, blkdiag(obj.Q_c, eye(obj.model.apx_nx)*obj.Q1), obj.R);
            K1 = K(1:obj.model.local_nx);
            K2 = K(obj.model.local_nx+1 : end);
            % initialize reward
            reward = 0;
            for k = 1 : obj.sim_N-1
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
        
        function [local_x_all, env_x_all, y_w_v_all, reward] = sim_original(obj, ini, noise_power, seed)
            if nargin < 4
                seed = rng();
            end
            if nargin < 3 || isempty(noise_power)
                noise_power = 1;
            end
            rng(seed);
            d_L = noise_power*randn(obj.sim_N, 2);
            local_x_all = zeros(obj.sim_N, obj.model.local_nx);
            env_x_all = zeros(obj.sim_N, obj.model.env_nx);
            y_w_v_all = zeros(obj.sim_N, obj.model.ny+obj.model.nw+obj.model.nv);
            % set initial
            local_x_all(1, :) = ini';
            % initialize reward
            reward = 0;
            for k = 1 : obj.sim_N-1
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

