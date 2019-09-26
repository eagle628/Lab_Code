classdef network_retro_by_actor_critic_with_eligibility_traces_episodic < RL_train
    % レトロフィット制御問題に対するRL_actor_criticのclass
    
    properties(Constant)
        Q1 = 1
        Q_c = diag([1,1])
        Q_r = diag([1,1])
        R = 1
        lambda_theta = 0.99
        lambda_omega = 0.99
        alpha = 0.05
        beta_mu = 0.01
        beta_sigma = 0.01
        alpha_f = 0.001 % fisher weghit %% invalid
        gamma = 0.99
        gamma2 = 0.9
        max_episode = 3.5e3
        snapshot = 100;
    end
    
    properties
        policy
        value
        sim_N
        t
        belief_N
    end
    
    methods
        function obj = network_retro_by_actor_critic_with_eligibility_traces_episodic(model, policy, value, Te, belief_N)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            obj.policy = policy;
            obj.value  =  value;
            obj.belief_N = belief_N;
        end
        
        function [local_x_all, mpc_u_all, rl_u_all, theta_mu_snapshot, theta_sigma_snapshot, w_snapshot, reward_history, varargout] = train(obj, ini, seed, varargin)
            if nargin < 3 || isempty(seed)
                seed = rng();
            end
            rng(seed)
            if nargout > 7
                varargout{1} = struct('cdata',[],'colormap',[]);
%                 varargout{2} = struct('cdata',[],'colormap',[]);
            end
            reward_history = zeros(obj.max_episode, 1);
%             w_history = zeros(obj.value.apx_function.N, obj.max_episode*obj.sim_N);
%             theta_mu_history = zeros(obj.policy.apx_function.N, obj.max_episode*obj.sim_N);
%             theta_sigma_history = zeros(obj.model.nu, obj.max_episode*obj.sim_N);
            % record point
            record_point = obj.snapshot:obj.snapshot:obj.max_episode;
            w_snapshot = zeros(obj.value.apx_function.N, length(record_point));
            theta_mu_snapshot = zeros(obj.policy.apx_function.N, length(record_point));
            theta_sigma_snapshot = zeros(obj.model.nu, length(record_point));
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
%             d_L = randn(obj.sim_N, 2);
            % params initialize
            % as much as possible (only information that we have)
%             w = set_w(obj, K);
            w = obj.value.get_params();
            theta_mu = obj.policy.get_params();
            theta_sigma = obj.policy.get_policy_sigma();
            % start episode learning
            record_idx = 1;
            belief_state = zeros(2, obj.belief_N*obj.model.ny);% belief state % 1 line: current state , 2 line: previous state
            for episode = 1 : obj.max_episode
                % memory reset
                local_x_all = nan(obj.sim_N, obj.model.local_nx);
                env_x_all = nan(obj.sim_N, obj.model.env_nx);
                rect_x_all = nan(obj.sim_N, obj.model.rect_nx);
                y_w_v_all = nan(obj.sim_N, obj.model.ny+obj.model.nw+obj.model.nv);
                rl_u_all = nan(obj.sim_N, obj.model.nu);
                mpc_u_all = nan(obj.sim_N, obj.model.nu);
                constraint_update = 0;
                obj.policy.initialize_memory();
                % episode initialize(eligibility)
                z_w = zeros(obj.value.apx_function.N, 1);
                z_theta_mu = zeros(obj.policy.apx_function.N, 1);
                z_theta_sigma = zeros(obj.model.nu, 1);
                Fix_Target_Network_params = w;
                zeta = 1;
                reward = 0;
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
                    belief_state(2, :) = belief_state(1, :);
                    belief_state(1, obj.model.ny+1:end) = belief_state(1, 1:end-obj.model.ny);
                    belief_state(1, 1:obj.model.ny) = rect_yw'; % momory store
                    % Get Reward r
%                     if abs(local_x_all(k,2) ) > 2
%                         r = -100;
%                     else
                        r = obj.reward([local_x_all(k+1, :), rect_x_all(k+1, :)], rl_u_all(k, :)+mpc_u_all(k, :));
% %                         r = 0;
%                     end
                    reward =  reward + obj.gamma^(k-1)*r;
                    % TD Erorr
                    tmp = strcmp(varargin, 'Fix-Target-Network');
                    if sum(tmp)
                        V_k1 = obj.value.est_value(belief_state(1, :), Fix_Target_Network_params);
                    else
                        V_k1 = obj.value.est_value(belief_state(1, :)); 
                    end
                    V_k0 = obj.value.est_value(belief_state(2, :));
                    delta = r + obj.gamma*V_k1 - V_k0;
                    tmp = strcmp(varargin, 'TD-Error-Clipping');
                    if sum(tmp)
                        if abs(delta) > varargin{find(tmp)+1}
                            delta = delta/abs(delta)*varargin{find(tmp)+1};
                        end
                    end
                    % eligibility traces update
                    z_w = obj.gamma*obj.lambda_theta*z_w + zeta*obj.value.value_grad(belief_state(2, :));
                    e_k1_mu = obj.policy.policy_grad_mu(rl_u_all(k, :), belief_state(2, :));
                    z_theta_mu = obj.gamma*obj.lambda_omega*z_theta_mu + zeta*e_k1_mu;
                    e_k1_sigma = obj.policy.policy_grad_sigma(rl_u_all(k, :), belief_state(2, :));
                    z_theta_sigma = obj.gamma*obj.lambda_omega*z_theta_sigma + zeta*e_k1_sigma;
                    % apx function update
                    if obj.policy.policy_constraint(theta_mu + obj.beta_mu*delta*z_theta_mu, obj.model, obj.belief_N)
                        w = w + obj.alpha*delta*z_w;
                        theta_mu = theta_mu + obj.beta_mu*delta*z_theta_mu;
                        theta_sigma = theta_sigma + obj.beta_sigma*delta*z_theta_sigma;
                        zeta = obj.gamma2*zeta;
                        obj.value.set_params(w);
                        obj.policy.set_params(theta_mu);
                        obj.policy.set_policy_sigma(theta_sigma);
                        constraint_update = constraint_update + 1;
                    end
                    % update Fix_Target_Network_params
                    tmp = strcmp(varargin, 'Fix-Target-Network');
                    if sum(tmp) && ~mod(episode, varargin{find(tmp)+1})
                        Fix_Target_Network_params = obj.value.get_params();
                    end
%                     if abs(local_x_all(k,2)) > 2% || abs(x_all(k,2)) > 4
%                         reward = -10;
%                         break;
%                     end
                end
                % record history
                if ~mod(episode, obj.snapshot)
                    w_snapshot(:, record_idx) = w; 
                    theta_mu_snapshot(:, record_idx) = theta_mu;
                    theta_sigma_snapshot(:, record_idx) = theta_sigma;
                    record_idx =  record_idx + 1;
                end
                reward_history(episode) = reward;
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
                if nargout > 7
                   varargout{1}(episode) = getframe(gcf);
               end
            end
        end
        
        function J = cost(obj, x_all, u_all)
           J = 0;
           for itr =  1 : obj.sim_N-1
               J = J + x_all(itr, :)*obj.Q_c*x_all(itr, :)' + u_all(itr, :)*obj.R*u_all(itr, :)';
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
            obj.policy.initialize_memory();
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

