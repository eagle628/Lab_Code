classdef REINFORCE_MC_episodic < RL_train
    % Episodic REINFORCE Monte-Carlo Policy Gradient episodic
    
    properties
        opt_policy
        gamma
        max_episode
        snapshot
        belief_sys
        render_enable
    end
    
    methods
        function obj = REINFORCE_MC_episodic(model, opt_policy, belief_sys)
            if nargin < 3 || isempty(belief_sys)
                belief_sys = observation_accumulater(model.ny, 1);
            end
            obj.model = model;
            obj.opt_policy = opt_policy;
            obj.gamma = 0.99;
            obj.max_episode = 10000;
            obj.snapshot = 1000;
            obj.belief_sys = belief_sys;
            obj.render_enable = true;
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
            for episode = 1 : obj.max_episode
                tic;
                % memory reset
                x_all = nan(obj.model.nx, sim_N);
                y_all = nan(obj.model.ny, sim_N);
                rl_u_all = nan(obj.model.nu, sim_N);
                rl_u_k_all = nan(obj.model.nu, sim_N);
                reward_all = nan(1, sim_N);
                % belief initialize
                belief_state_all = nan(size(obj.belief_sys.A, 1), sim_N);
                % set initial
                obj.model.initialize(ini_set(:, episode));
                 % episode initialize
                obj.opt_policy.initialize(episode);
                obj.belief_sys.initialize();
                reward = 0;
                % belief initialize
                y_all(:, 1) = obj.model.observe();
                x_all(:, 1) = obj.model.state;
                belief_state_all(:, 1) =  obj.belief_sys.estimate(y_all(:, 1));
                for k = 1 : sim_N-1
                    [rl_u_all(:, k), rl_u_k_all(:,k)] = obj.opt_policy.target.predict(belief_state_all(:, k), true);
                    % next step
                    [y_all(:, k+1), r] = obj.model.dynamics(rl_u_all(:, k));
                    x_all(:, k+1) = obj.model.state;
                    reward_all(k) = r;
                    reward =  reward + obj.gamma^(k-1)*r;
                    % belief upadate
%                     obj.belief_sys.update_internal_state(y_all(:, k), rl_u_all(:, k));% state estimator  internal state update
                    obj.belief_sys.update_internal_state(y_all(:, k), rl_u_k_all(:, k));% state estimator  internal state update
                    belief_state_all(:, k+1) = obj.belief_sys.estimate(y_all(:, k+1));
                end
                % parameter update
                for itr1 = obj.belief_sys.accumulate_N : sim_N-1
                    delta = 0;
                    for itr2 = itr1 : sim_N
                        delta = delta + obj.gamma^(itr2-itr1)*reward_all(itr2);
                    end
                    data.delta = delta;
                    data.state = belief_state_all(:, itr1);
                    data.pre_input = rl_u_all(:, k);
                    data.pre_input_mu = rl_u_k_all(:, itr1);
                    % opt
                    obj.opt_policy.opt(data);
                end
                % record history
                if ~mod(episode, obj.snapshot) 
                    policy_snapshot{record_idx} = obj.opt_policy.target.get_params();
                    record_idx =  record_idx + 1;
                end
                history.reward(episode) = reward;
                history.policy_counter(episode) = obj.opt_policy.counter;
                if obj.render_enable
                    obj.render(t, x_all, y_all, history.reward, episode);
                end
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
        
        function render(obj, t, x_all, y_all, reward_history, episode, varargin)
            subplot(2,1,1)
            plot(t, y_all)
            title(['\fontsize{16}','Episode-',num2str(episode)])
            ylabel('y')
            grid on
            lim = max(max(abs(y_all),[],'omitnan'));
            if isempty(lim) || lim == 0
                lim = 1;
            end
            ylim([-lim, lim]);
            subplot(2,1,2)
            plot(nonzeros(reward_history),'-b')
            ylabel('Culumative Reward')
            drawnow
        end
    end
end

