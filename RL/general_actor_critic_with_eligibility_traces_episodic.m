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
        max_episode = 1e4
        snapshot = 100
    end
    
    properties
        policy
        value
        sim_N
        t
    end
    
    methods
        function obj = general_actor_critic_with_eligibility_traces_episodic(model, policy, value, Te)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            obj.policy = policy;
            obj.value  =  value;
        end
        
        function [x_all, mpc_u_all, rl_u_all, theta_mu_snapshot, theta_sigma_snapshot, w_snapshot, reward_history, varargout] = train(obj, ini, seed, varargin)
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
            theta_sigma_snapshot = zeros(obj.model.nu, length(record_point));
            % calculate MBC gain
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            tmp = strcmp(varargin, 'parallel');
            if sum(tmp)
                if  strcmp(varargin{find(tmp)+1},'off')
                    K = zeros(size(K));
                end
            end
            % params initialize
            % as much as possible (only information that we have)
%             w = set_w(obj, K);
            w = obj.value.get_params();
            theta_mu = obj.policy.get_params();
            theta_sigma = obj.policy.get_policy_sigma();
            % start episode learning
            record_idx = 1;
            for episode = 1 : obj.max_episode
                % episode initialize
                z_w = zeros(obj.value.apx_function.N, 1);
                z_theta_mu = zeros(obj.policy.apx_function.N, 1);
                z_theta_sigma = zeros(obj.model.nu, 1);
                zeta = 1;
                reward = 0;
                % memory reset
                x_all = nan(obj.sim_N, obj.model.true_nx);
                y_all = nan(obj.sim_N, obj.model.ny);
                rl_u_all = nan(obj.sim_N, obj.model.nu);
                mpc_u_all = nan(obj.sim_N, obj.model.nu);
                % set initial
                x_all(1, :) = [rand(1) - 0.5, 0];
%                 x_all(1, :) = ini';
                % belief initialize
                for k = 1 : obj.sim_N-1
                    % MBC input
                    mpc_u_all(k, :) = -K*x_all(k, :)';
%                    mpc_u_all(k, :) = 0;
                    % RL input
                    tmp = strcmp(varargin, 'Input-Clipping');
                    if sum(tmp)
                        rl_u_all(k, :) = obj.policy.stocastic_policy(x_all(k, :), theta_mu, 'Input-Clipping', varargin{find(tmp)+1});
                    else
                        rl_u_all(k, :) = obj.policy.stocastic_policy(x_all(k, :), theta_mu);
                    end
                    %  observe S_(k+1)
                   [ne_x, y] = obj.model.dynamics(x_all(k, :)', rl_u_all(k, :) + mpc_u_all(k, :));
                    x_all(k+1, :) = ne_x';
                    y_all(k, :) = y';
                    % Get Reward r
%                     if abs(x_all(k,1)) > 0.5 || abs(x_all(k,2)) > 4
%                         r = -3;
%                     else
                        r = obj.reward(x_all(k+1, :), rl_u_all(k, :)+mpc_u_all(k, :));
%                     end
                    reward =  reward + obj.gamma^(k-1)*r;
                    % TD Erorr
                    V_k1 = obj.value.est_value(x_all(k+1, :), w);
                    V_k0 = obj.value.est_value(x_all(k, :), w);
                    delta = r + obj.gamma*V_k1 - V_k0;
                    tmp = strcmp(varargin, 'TD-Error-Clipping');
                    if sum(tmp)
                        delta = -delta/delta*varargin{find(tmp)+1};
                    end
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
%                     figure(3)
%                     stem(w)
%                     drawnow
                    if abs(x_all(k,1)) > 0.5
                        break;
                    end
                end
%                close 3
                % record history
                if ~mod(episode, obj.snapshot) 
                    w_snapshot(:, record_idx) = w; 
                    theta_mu_snapshot(:, record_idx) = theta_mu;
                    theta_sigma_snapshot(:, record_idx) = theta_sigma;
                    record_idx =  record_idx + 1;
                end
                reward_history(episode) = reward;
                cost_history(episode) = obj.cost(x_all, rl_u_all+mpc_u_all);
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
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            tmp = strcmp(varargin, 'parallel');
            if sum(tmp)
                if  strcmp(varargin{find(tmp)+1},'off')
                    K = zeros(size(K));
                end
            end
            for itr = 1 : obj.sim_N-1
%                 u_rl = obj.policy.determistic_policy(x_all(itr, :), theta);
                u_rl = obj.policy.stocastic_policy(x_all(itr, :), theta);
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

