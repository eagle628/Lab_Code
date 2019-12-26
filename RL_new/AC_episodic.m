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
        render_enable
        value_pretraining_period
    end
    
    methods
        function obj = AC_episodic(model, opt_policy, opt_value, belief_sys)
            if nargin < 4 || isempty(belief_sys)
                belief_sys = observation_accumulater(model.ny, 1);
            end
            obj.model = model;
            obj.opt_policy = opt_policy;
            obj.opt_value  = opt_value;
            obj.gamma = 0.99;
            obj.max_episode = 10000;
            obj.snapshot = 1000;
            obj.belief_sys = belief_sys;
            obj.fixed_apx_function_period = 100;
            obj.render_enable = true;
            obj.value_pretraining_period = 0;
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
            history.reward = zeros(obj.max_episode, 1);
            history.policy_counter = zeros(obj.max_episode, 1);
            history.value_counter = zeros(obj.max_episode, 1);
%             history.delta = cell(obj.max_episode, 1);
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
            belief_state = zeros(size(obj.belief_sys.A, 1), 2);% belief state % 1 line: current state , 2 line: previous state
            for episode = 1 : obj.max_episode
                tic;
                % memory reset
                x_all = nan(obj.model.nx, sim_N);
                y_all = nan(obj.model.ny, sim_N);
                rl_u_all = nan(obj.model.nu, sim_N);
%                 delta_history = cell(sim_N-1, 1);
                % belief initialize
                belief_state = zeros(size(belief_state));
                % set initial
                obj.model.initialize(ini_set(:, episode));
                 % episode initialize
                obj.opt_policy.initialize(episode);
                obj.opt_value.initialize(episode);
                obj.belief_sys.initialize();
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
                belief_state(:, 1) =  obj.belief_sys.estimate(y_all(:, 1));
                % 
                grad = 0;
                x_ = obj.opt_policy.target.get_params();
                x_ = x_(1:end-1);
                input_ = zeros(size(rl_u_all));
                delta_set = zeros(sim_N,1);
%                 rng(6)
                for k = 1 : sim_N-1
                    actor_belief = zeros(size(belief_state(:, 1)));
                    actor_belief(1:obj.model.ny) = belief_state(1:obj.model.ny, 1);
                    [rl_u_all(:, k), rl_u_k] = obj.opt_policy.target.predict(belief_state(:, 1), true);
%                     [rl_u_all(:, k), rl_u_k] = obj.opt_policy.target.predict(actor_belief, true);
                    input_(:, k) = rl_u_all(:, k) - rl_u_k;
                    % next step
                    [y_all(:, k+1), r] = obj.model.dynamics(rl_u_all(:, k));
                    x_all(:, k+1) = obj.model.state;
                    reward =  reward + obj.gamma^(k-1)*r;
                    % belief upadate
                    belief_state(:, 2) = belief_state(:, 1);% save previous belief state
%                     obj.belief_sys.update_internal_state(y_all(:, k), rl_u_all(:, k));% state estimator  internal state update
                    obj.belief_sys.update_internal_state(y_all(:, k), rl_u_k);% state estimator  internal state update
                    belief_state(:, 1) = obj.belief_sys.estimate(y_all(:, k+1));
                    % TD Erorr
                    V_k1 = obj.opt_value.target.predict(belief_state(:, 1), false);
                    V_k0 = obj.opt_value.target.predict(belief_state(:, 2), true);
                    delta = r + obj.gamma*V_k1 - V_k0;
                    % parameter update
                    data.delta = delta;
                    data.state = belief_state(:, 2);
%                     data.state = actor_belief;
                    data.pre_input = rl_u_all(:, k);
                    data.pre_input_mu = rl_u_k;
                    if isa(data.delta, 'double')
                        delta_set(k) = data.delta;
                    else
                        delta_set(k) = double(py.numpy.squeeze(data.delta.data).tolist);
                    end
                    if k > obj.belief_sys.accumulate_N
                        if episode > obj.value_pretraining_period
                            grad = grad + obj.opt_policy.opt(data);
                        end
%                     data.state = belief_state(:, 2);
%                         if episode < 1
%                             obj.opt_value.opt(data);
%                         end
                    end
                    % for pendulum
%                     if abs(x_all(1, k)) > 0.5
%                         reward = -30;
%                         break;
%                     end
                    % save delta
%                     delta_history{k} = data.delta;
                end
                counter_ = obj.opt_policy.counter;
%                 save(sprintf('test_rl_grad_n%d_belief_%d_epi%d',obj.model.net.N,obj.belief_sys.accumulate_N, episode), 'x_', 'grad' ,'input_', 'counter_');
                % record history
                if ~mod(episode, obj.snapshot) 
                    value_snapshot{record_idx}  = obj.opt_value.target.get_params();
                    policy_snapshot{record_idx} = obj.opt_policy.target.get_params();
                    record_idx =  record_idx + 1;
                end
                history.reward(episode) = reward;
                history.policy_counter(episode) = obj.opt_policy.counter;
                history.value_counter(episode) = obj.opt_value.counter;
%                 history.delta{episode} = delta_history;
                if obj.render_enable
                    obj.render(t, x_all, y_all, history.reward, episode, delta_set);
                end
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
            reward = 0;
            obj.opt_policy.initialize(0);
            obj.belief_sys.initialize();
            belief_state =  obj.belief_sys.estimate(y_all(:, 1));
            for k = 1 : sim_N-1
                rl_u_all(:, k) = obj.opt_policy.target.predict(belief_state, false);
                [y_all(:, k+1), r] = obj.model.dynamics(rl_u_all(:, k));
                x_all(:, k+1) = obj.model.state;
                obj.belief_sys.update_internal_state(y_all(:, k), rl_u_all(:, k));
                belief_state =  obj.belief_sys.estimate(y_all(:, k+1));
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
            obj.opt_policy.initialize(0);
            for k = 1 : sim_N-1
                mpc_u_all(:, k) = -K*x_all(:, k);
                [y_all(:, k+1), r] = obj.model.dynamics(mpc_u_all(:, k));
                x_all(:, k+1) = obj.model.state;
                reward = reward + obj.gamma^(k-1)*r;
            end
        end
        
        function render(obj, t, x_all, y_all, reward_history, episode, delta)
            subplot(3,2,1)
            plot(t, y_all)
            title(['\fontsize{16}','Episode-',num2str(episode)])
            ylabel('y')
            grid on
            lim = max(max(abs(y_all),[],'omitnan'));
            if isempty(lim) || lim == 0
                lim = 1;
            end
            ylim([-lim, lim]);
            subplot(3,2,3)
            plot(nonzeros(reward_history),'-b')
            ylabel('Culumative Reward')
            subplot(3,2,5)
            plot(episode, norm(delta),'r*')
            ylabel('TD-error')
            hold on
            drawnow
            % %
            subplot(3,2,[2,4,6])
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

