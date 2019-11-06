classdef AC_episodic < RL_train
    % Episodic Actor Critic trainer
    
    properties
        opt_policy
        opt_value
        gamma
        max_episode
        snapshot
        belief_sys
        fixed_apx_function_period
    end
    
    methods
        function obj = AC_episodic(model, opt_policy, opt_value, belief_sys)
            if nargin < 4
                belief_sys = [zeros(model.ny), eye(model.ny)];
            end
            obj.model = model;
            obj.opt_policy = opt_policy;
            obj.opt_value  = opt_value;
            obj.gamma = 0.99;
            obj.max_episode = 10000;
            obj.snapshot = 1000;
            obj.belief_sys = belief_sys;
            obj.fixed_apx_function_period = 100;
        end
        
        function [x_all, rl_u_all, policy_snapshot, value_snapshot, history] = train(obj, ini_set, Te, seed, varargin)
            sim_N = Te/obj.model.Ts + 1;
            t = (0:obj.model.Ts:Te)';
            if nargin < 3 || isempty(seed)
                seed = rng();
            end
            rng(seed)
            % history
            history = struct();
            history.reward_= zeros(obj.max_episode, 1);
            history.policy_conter = zeros(obj.max_episode, 1);
            history.value_conter = zeros(obj.max_episode, 1);
            history.delta = cell(obj.max_episode, 1);
            % record point
            record_point = obj.snapshot:obj.snapshot:obj.max_episode;
            value_snapshot = cell(1, length(record_point));
            policy_snapshot = cell(1, length(record_point));
            % 
            data = struct();
            data.model = obj.model;
            data.belief_sys = obj.belief_sys;
            data.gamma = obj.gamma;
            data.state = [];
            data.delta = [];
            data.pre_input = [];
            data.pre_input_mu = [];
            % start episode learning
            record_idx = 1;
            belief_state = zeros(size(obj.belief_sys, 1), 2);% belief state % 1 line: current state , 2 line: previous state
            for episode = 1 : obj.max_episode
                tic;
                % memory reset
                x_all = nan(obj.model.nx, sim_N);
                y_all = nan(obj.model.ny, sim_N);
                rl_u_all = nan(obj.model.nu, sim_N);
                % belief initialize
                belief_state = zeros(size(belief_state));
                % set initial
                obj.model.initialize(ini_set(:, episode));
                 % episode initialize
                obj.opt_policy.initialize(episode);
                obj.opt_value.initialize(episode);
                reward = 0;
                % fixed apx function
                if isprop(obj.opt_value.target,'fixed_apx_function_enable') && obj.opt_value.target.fixed_apx_function_enable
                    if ~mod(episode, obj.fixed_apx_function_period)
                            obj.opt_value.target.fixed_apx_function_update();
                    end
                end
                % belief initialize
                y_all(:, 1) = obj.model.observe();
                x_all(:, 1) = obj.model.state;
                belief_state(:, 1) =  obj.belief_sys*[belief_state(:, 1); y_all(:, 1)];
                for k = 1 : sim_N-1
                    [rl_u_all(:, k), rl_u_k] = obj.opt_policy.target.predict(belief_state(:, 1), true);
                    % next step
                    [y_all(:, k+1), r] = obj.model.dynamics(rl_u_all(:, k));
                    x_all(:, k+1) = obj.model.state;
                    reward =  reward + obj.gamma^(k-1)*r;
                    % belief upadate
                    belief_state(:, 2) = belief_state(:, 1);
                    belief_state(:, 1) = obj.belief_sys*[belief_state(:, 1); y_all(:, k+1)]; % momory store
                    % TD Erorr
                    V_k1 = obj.opt_value.target.predict(belief_state(:, 1), false);
                    V_k0 = obj.opt_value.target.predict(belief_state(:, 2), true);
                    delta = r + obj.gamma*V_k1 - V_k0;
                    % parameter update
                    data.delta = delta;
                    data.state = belief_state(:, 2);
                    data.pre_input = rl_u_all(:, k);
                    data.pre_input_mu = rl_u_k;
                    obj.opt_policy.opt(data);
                    obj.opt_value.opt(data);
%                     if abs(x_all(1, k)) > 0.5
%                         reward = -30;
%                         break;
%                     end
                end
                % record history
                if ~mod(episode, obj.snapshot) 
                    value_snapshot{record_idx}  = obj.opt_value.target.get_params();
                    policy_snapshot{record_idx} = obj.opt_policy.target.get_params();
                    record_idx =  record_idx + 1;
                end
                history.reward(episode) = reward;
                history.policy_counter(episode) = obj.opt_policy.counter;
                history.value_counter(episode) = obj.opt_value.counter;
                history.delta{episode} = data.delta;
                obj.render(t, x_all, y_all, history.reward, episode);
                disp(strcat('Episode-',num2str(episode),' : value  constraint update times : ', num2str(obj.opt_value.counter) ,'/',num2str(k)))
                disp(strcat('Episode-',num2str(episode),' : policy constraint update times : ', num2str(obj.opt_policy.counter) ,'/',num2str(k)))
                timer = toc;
                fprintf('This epoch %f[s], Estimated time to finish:%f [h].\n',timer, timer*(obj.max_episode-episode)/3600)
            end
        end
        
        function [x_all, y_all, rl_u_all, t, reward] = sim(obj, ini, Te, varargin)
            sim_N = Te/obj.model.Ts + 1;
            t = (0:obj.model.Ts:Te)';
            x_all = nan(obj.model.nx, sim_N);
            y_all = nan(obj.model.ny, sim_N);
            rl_u_all = nan(obj.model.nu, sim_N);
            % set initial
            obj.model.initialize(ini);
            y_all(:, 1) = obj.model.observe();
            x_all(:, 1) = obj.model.state;
            belief_state = zeros(size(obj.belief_sys, 1), 1);
            belief_state =  obj.belief_sys*[belief_state; y_all(:, 1)];
            reward = 0;
            obj.opt_policy.initialize();
            for k = 1 : sim_N-1
                rl_u_all(:, k) = obj.opt_policy.target.predict(belief_state, false);
                [y_all(:, k+1), r] = obj.model.dynamics(rl_u_all(:, k));
                x_all(:, k+1) = obj.model.state;
                belief_state =  obj.belief_sys*[belief_state; y_all(:, k+1)];
                reward = reward + obj.gamma^(k-1)*r;
            end
        end
        
        function [x_all, y_all, mpc_u_all, t, reward] = sim_lqrcontroller(obj, ini, Te, varargin)
            sim_N = Te/obj.model.Ts + 1;
            t = (0:obj.model.Ts:Te)';
            K = dlqr(obj.model.A, obj.model.B, obj.model.Q, obj.model.R);
            x_all = nan(obj.model.nx, sim_N);
            y_all = nan(obj.model.ny, sim_N);
            mpc_u_all = nan(obj.model.nu, sim_N);
            % set initial
            obj.model.initialize(ini);
            y_all(:, 1) = obj.model.observe();
            x_all(:, 1) = obj.model.state;
            reward = 0;
            obj.opt_policy.initialize();
            for k = 1 : sim_N-1
                mpc_u_all(:, k) = -K*x_all(:, k);
                [y_all(:, k+1), r] = obj.model.dynamics(mpc_u_all(:, k));
                x_all(:, k+1) = obj.model.state;
                reward = reward + obj.gamma^(k-1)*r;
            end
        end
        
        function render(obj, t, x_all, y_all, reward_history, episode, update_chance)
            subplot(2,2,1)
            plot(t, y_all)
            title(['\fontsize{16}','Episode-',num2str(episode)])
            ylabel('y')
            grid on
            lim = max(max(abs(y_all),[],'omitnan'));
            if isempty(lim) || lim == 0
                lim = 1;
            end
            ylim([-lim, lim]);
            subplot(2,2,3)
            plot(nonzeros(reward_history),'-b')
            ylabel('Culumative Reward')
            drawnow
            disp(strcat('Episode-',num2str(episode),' : value  constraint update times : ', num2str(obj.opt_value.counter) ,'/',num2str(update_chance)))
            disp(strcat('Episode-',num2str(episode),' : policy constraint update times : ', num2str(obj.opt_policy.counter) ,'/',num2str(update_chance)))
            timer = toc;
            fprintf('This epoch %f[s], Estimated time to finish:%f [h].\n',timer, timer*(obj.max_episode-episode)/3600)
            % %
            subplot(2,2,[2,4])
            [X,Y] = meshgrid(-.5:.1:.5, -2:.4:2);
            mesh_size = size(X, 1);
            XY = zeros(2*mesh_size, mesh_size);
            XY(1:2:end, :) = X;
            XY(2:2:end, :) = Y;
            XY = mat2cell(XY, 2*ones(1,mesh_size), ones(1,mesh_size));
            Z = cellfun(@(x)obj.opt_value.target.predict(x), XY);
            mesh(X,Y,Z)
            xlabel('x_1')
            ylabel('x_2')
            zlabel('Value')
            title(['\fontsize{16}','Episode-',num2str(episode)]);
        end
        
% %         function w = set_w(obj, K, w0)
% %             if nargin < 3
% %                 w0 = randn(obj.value.apx_function.N, 1);
% %             end
% %             [X,Y] = meshgrid(-2:.1:2, -2:.1:2);
% %             mesh_size = size(X, 1);
% %             Z = zeros(mesh_size, mesh_size);
% %             for itr1 = 1 : mesh_size
% %                for itr2 = 1 : mesh_size
% %                     Z(itr1, itr2) = obj.cumulative_reward(K, [X(itr1,itr2);Y(itr1,itr2)]);
% %                end
% %             end
% %             state = [vec(X), vec(Y)];
% % %             Z = -1*ones(size(Z));
% %             target = vec(Z);
% % 
% %             options = optimoptions('lsqnonlin','Display','iter','SpecifyObjectiveGradient',true);
% %             w = lsqnonlin(@(w)obj.apx_cost_function(state, target, w),w0,[],[],options);
% % % %             
% % % %             [X,Y] = meshgrid(-0.5:.1:0.5, -2:.4:2);
% % % %             mesh_size = size(X, 1);
% % % %             Z = zeros(mesh_size, mesh_size);
% % % %             for itr1 = 1 : mesh_size
% % % %                for itr2 = 1 :mesh_size
% % % %                     Z(itr1, itr2) = obj.value.est_value([X(itr1,itr2),Y(itr1,itr2)], w);
% % % %                end
% % % %             end
% % % %             mesh(X,Y,Z)
% %             
% %         end
% %         
% %         function [r ,apx_all, u_all] = cumulative_reward(obj, K, ini)
% %             apx_all = zeros(obj.sim_N, obj.model.apx_nx);
% %             y_all = zeros(obj.sim_N, obj.model.ny);
% %             u_all = zeros(obj.sim_N, obj.model.nu);
% %             r = 0;
% %             % set initial
% %             apx_all(1, :) = ini';
% %             for k = 1 : obj.sim_N-1
% %                 u_all(k, :) = (-K*apx_all(k, :)')';
% %                 [ne_x, y] = obj.model.approximate_dynamics(apx_all(k, :)', u_all(k, :)');
% %                 apx_all(k+1, :) = ne_x';
% %                 y_all(k, :) = y';
% %                 r = r + obj.gamma^(k-1)*obj.reward(apx_all(k+1, :), u_all(k, :));
% %             end
% %         end
% %         
% %         function [error, grad] = apx_cost_function(obj, x, target, w)
% %             apx = zeros(size(target));
% %             for itr = 1 : length(target)
% %                 apx(itr) = obj.value.est_value(x(itr, :), w);
% %             end
% %             error = apx - target;
% %             if nargout > 1
% %                 grad = zeros(length(target), obj.value.apx_function.N);
% %                 for itr = 1 : length(target)
% %                     grad(itr, :) = obj.value.value_grad(x(itr, :));
% %                 end
% %             end
% %         end
    end
end

