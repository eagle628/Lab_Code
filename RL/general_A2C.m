classdef general_A2C < RL_train
    % 一般的な制御問題に対するRL_actor_criticのclass
    
    properties(Constant)
        Q = diag([10,1])
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
    end
    
    methods
        function obj = general_A2C(model, policy, value, Te, advantage_N)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            obj.policy = policy;
            obj.value  =  value;
            if nargin < 5
                advantage_N = 5;
            end
            obj.advantage_N = advantage_N;
        end
        
        function [x_all, rl_u_all, theta_mu_snapshot, w_snapshot, reward_history, varargout] = train(obj, ini, seed, varargin)
            if nargin < 3 || isempty(seed)
                seed = rng();
            end
            rng(seed)
            rng('shuffle')
            if nargout > 7
                varargout{1} = struct('cdata',[],'colormap',[]);
                varargout{2} = struct('cdata',[],'colormap',[]);
            end
            cost_history = zeros(obj.max_episode, 1);
            reward_history = zeros(obj.max_episode, 1);
            % record point
            record_point = obj.snapshot:obj.snapshot:obj.max_episode;
            w_snapshot = zeros(obj.value.apx_function.N, length(record_point));
            theta_mu_snapshot = zeros(obj.policy.apx_function.N, length(record_point));
            % params initialize
            w = obj.value.get_params();
            theta_mu = obj.policy.get_params();
            % start episode learning
            record_idx = 1;
            for episode = 1 : obj.max_episode
                reward = 0;
                % memory reset
                x_all = nan(obj.sim_N, obj.model.true_nx);
                y_all = nan(obj.sim_N, obj.model.ny);
                rl_u_all = nan(obj.sim_N, obj.model.nu);
                % set initial
                x_all(1, :) = [rand(1) - 0.5, 0];
%                 x_all(1, :) = ini';
                % 
                for k = 1 : obj.sim_N-1
                    % RL input
                    tmp = strcmp(varargin, 'Input-Clipping');
                    if sum(tmp)
                        rl_u_all(k, :) = obj.policy.stocastic_policy(x_all(k, :), theta_mu, 'Input-Clipping', varargin{find(tmp)+1});
                    else
                        rl_u_all(k, :) = obj.policy.stocastic_policy(x_all(k, :), theta_mu);
                    end
                    %  observe S_(k+1)
                   [ne_x, y] = obj.model.dynamics(x_all(k, :)', rl_u_all(k, :));
                    x_all(k+1, :) = ne_x';
                    y_all(k, :) = y';
                    reward =  reward + obj.gamma^(k-1)*obj.reward(x_all(k+1, :), rl_u_all(k, :));
                    if ~mod(k, obj.advantage_N)
                        % reset accumulate gradients
                        d_theta_mu = zeros(size(theta_mu));
                        d_w        = zeros(size(w));
                        backward_reward = obj.value.est_value(x_all(k+1, :), w);
                        for iter1 = k : -1 : k - obj.advantage_N+1
                            % Get Reward r
                            r = obj.reward(x_all(iter1+1, :), rl_u_all(iter1, :));
                            backward_reward = r + obj.gamma*backward_reward;
                            advantage = backward_reward - obj.value.est_value(x_all(iter1, :)); 
                            % accumulates gradients
                            d_theta_mu = d_theta_mu + obj.policy.policy_grad_mu(rl_u_all(iter1, :), x_all(iter1, :))*advantage;
                            d_w        = d_w        + obj.value.value_grad(x_all(iter1, :))*advantage;
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
                        end
                        % update params
                        obj.value.set_params(w+obj.alpha*d_w);
                        obj.policy.set_params(theta_mu+obj.beta*d_theta_mu);
                        w = obj.value.get_params();
                        theta_mu = obj.policy.get_params();
                    end
                    if abs(x_all(k,1)) > 0.5
                        break;
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
                cost_history(episode) = obj.cost(x_all, rl_u_all);
                figure(1)
                callback_RL(episode, obj.t, x_all, cost_history, reward_history)
                if nargout > 7
                   varargout{1}(episode) = getframe(gcf);
                end
                figure(2)
                [X,Y] = meshgrid(-0.5:.1:0.5, -2:.4:2);
                mesh_size = size(X, 1);
                XY = zeros(mesh_size, mesh_size*2);
                XY(:, 1:2:end) = X;
                XY(:,2:2:end) = Y;
                XY = mat2cell(XY, ones(1,mesh_size), 2*ones(1,mesh_size));
                Z = cellfun(@(x)obj.value.est_value(x), XY);
                mesh(X,Y,Z)
                xlabel('x_1')
                ylabel('x_2')
                zlabel('Value')
                title(strcat('episode-',num2str(episode)));
                if nargout > 7
                   varargout{2}(episode) = getframe(gcf);
                end
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
        
        function x_all = sim(obj, ini, theta, varargin)
            if nargin < 3 || isempty(theta)
                theta = obj.policy.get_params();
            end
            x_all = zeros(obj.sim_N, obj.model.true_nx);
            x_all(1, :) = ini';
            for itr = 1 : obj.sim_N-1
%                 u_rl = obj.policy.determistic_policy(x_all(itr, :), theta);
                u_rl = obj.policy.stocastic_policy(x_all(itr, :), theta);
                ne_x = (obj.model.dynamics(x_all(itr,:)', u_rl))';
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

