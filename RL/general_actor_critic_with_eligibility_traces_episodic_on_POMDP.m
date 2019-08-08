classdef general_actor_critic_with_eligibility_traces_episodic_on_POMDP < RL_train
    % 一般的な制御問題に対するRL_actor_criticのclass
    
    properties(Constant)
        lambda_theta = 0.99
        lambda_omega = 0.99
        alpha = 0.0005
        beta_mu = 0.0001
        beta_sigma = 0.01
        gamma = 0.99
        gamma2 = 0.9
        max_episode = 4e4
        snapshot = 100
    end
    
    properties
        policy
        value
        sim_N
        t
        belief_N
        Q
        R
    end
    
    methods
        function obj = general_actor_critic_with_eligibility_traces_episodic_on_POMDP(model, policy, value, Te, belief_N)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            obj.policy = policy;
            obj.value  =  value;
            obj.belief_N = belief_N;
            obj.Q = blkdiag(1,zeros(belief_N*model.ny-1));
            obj.R = 0;
        end
        
        function [x_all, rl_u_all, theta_mu_snapshot, theta_sigma_snapshot, w_snapshot, reward_history, varargout] = train(obj, ini, seed, varargin)
            if nargin < 3 || isempty(seed)
                seed = rng();
            end
            rng(seed)
            if nargout > 6
                varargout{1} = struct('cdata',[],'colormap',[]);
%                 varargout{2} = struct('cdata',[],'colormap',[]);
            end
            cost_history = zeros(obj.max_episode, 1);
            reward_history = zeros(obj.max_episode, 1);
            % record point
            record_point = obj.snapshot:obj.snapshot:obj.max_episode;
            w_snapshot = zeros(obj.value.apx_function.N, length(record_point));
            theta_mu_snapshot = zeros(obj.policy.apx_function.N, length(record_point));
            theta_sigma_snapshot = zeros(obj.model.nu, length(record_point));
            % params initialize
            % as much as possible (only information that we have)
            w = obj.value.get_params();
            theta_mu = obj.policy.get_params();
            theta_sigma = obj.policy.get_policy_sigma();
            Fix_Target_Network_params = w;
            % start episode learning
            record_idx = 1;
            % observe belief state
            belief_state = zeros(2, obj.belief_N*obj.model.ny);% belief state % 1 line: current state , 2 line: previous state
            for episode = 1 : obj.max_episode
                % memory reset
                x_all = nan(obj.sim_N, obj.model.true_nx);
                y_all = nan(obj.sim_N, obj.model.ny);
                rl_u_all = nan(obj.sim_N, obj.model.nu);
                % episode initialize
                z_w = zeros(obj.value.apx_function.N, 1);
                z_theta_mu = zeros(obj.policy.apx_function.N, 1);
                z_theta_sigma = zeros(obj.model.nu, 1);
                zeta = 1;
                reward = 0;
                % set episode initial
%                 x_all(1, :) = ini';
                x_all(1, :) = [rand(1) - 0.5, 0];
                % belief initialize
                belief_state = zeros(size(belief_state));
                for k = 1 : obj.sim_N-1
                    % RL input
                    tmp = strcmp(varargin, 'Input-Clipping');
                    if sum(tmp)
                        rl_u_all(k, :) = obj.policy.stocastic_policy(belief_state(1, :), [], 'Input-Clipping', varargin{find(tmp)+1});
                    else
                        rl_u_all(k, :) = obj.policy.stocastic_policy(belief_state(1, :), []);
                    end
                    %  observe S_(k+1)
                   [ne_x, y] = obj.model.dynamics(x_all(k, :)', rl_u_all(k, :));
                    x_all(k+1, :) = ne_x';
                    y_all(k, :) = y';
                    % belief upadate
                    belief_state(2, :) = belief_state(1, :);
                    belief_state(1, obj.model.ny+1:end) = belief_state(1, 1:end-obj.model.ny);
                    belief_state(1, 1:obj.model.ny) = y'; % momory store
                    % Get Reward r
%                     if abs(x_all(k,1)) > 0.5 || abs(x_all(k,2)) > 4
%                         r = -3;
%                     else
                        r = obj.reward(belief_state(1, :), rl_u_all(k, :));
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
                    w = w + obj.alpha*delta*z_w;
                    theta_mu = theta_mu + obj.beta_mu*delta*z_theta_mu;
                    theta_sigma = theta_sigma + obj.beta_sigma*delta*z_theta_sigma;
                    zeta = obj.gamma2*zeta;
                    obj.value.set_params(w);
                    obj.policy.set_params(theta_mu);
                    obj.policy.set_policy_sigma(theta_sigma);
                    % update Fix_Target_Network_params
                    tmp = strcmp(varargin, 'Fix-Target-Network');
                    if sum(tmp) && ~mod(episode, varargin{find(tmp)+1})
                        Fix_Target_Network_params = obj.value.get_params();
                    end
%                     figure(3)
%                     stem(w)
%                     drawnow
                    if abs(x_all(k,1)) > 0.5% || abs(x_all(k,2)) > 4
                        reward = -10;
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
%                cost_history(episode) = obj.cost(x_all, rl_u_all);% not calculation
               figure(1)
               callback_RL(episode, obj.t, x_all, cost_history, reward_history)
               if nargout > 6
                   varargout{1}(episode) = getframe(gcf);
               end
%                figure(2)
%                [X,Y] = meshgrid(-0.5:.1:0.5, -0.5:.1:0.5);
%                mesh_size = size(X, 1);
%                XY = zeros(mesh_size, mesh_size*2);
%                XY(:,1:2:end) = X;
%                XY(:,2:2:end) = Y;
%                XY = mat2cell(XY, ones(1,mesh_size), 2*ones(1,mesh_size));
%                Z = cellfun(@(x)obj.value.est_value(x), XY);
%                mesh(X,Y,Z)
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
        
        function x_all = sim(obj, ini, theta)
            if nargin < 3 || isempty(theta)
                theta = obj.policy.get_params();
            end
            x_all = zeros(obj.sim_N, obj.model.true_nx);
            x_all(1, :) = ini';
            % observe belief state
            belief_state = zeros(1, obj.belief_N);
            for itr = 1 : obj.sim_N-1
                u_rl = obj.policy.determistic_policy(belief_state, theta);
%                 u_rl = obj.policy.stocastic_policy(belief_state, theta);
                [ne_x, y] = obj.model.dynamics(x_all(itr,:)', u_rl);
                x_all(itr+1, :) = ne_x';
                belief_state(1, 1+1:end) = belief_state(1, 1:end-1);
                belief_state(1, 1) = y'; % momory store
            end
        end
        
%         function [x_all] = sim_lqrcontroller(obj, ini)
%             x_all = zeros(obj.sim_N, 2);
%             x_all(1, :) = ini';
%             K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
%             for itr = 1 : obj.sim_N-1
%                 u_mbc = -K*x_all(itr, :)';
%                 x_all(itr+1, :) = (obj.model.dynamics(x_all(itr,:)', u_mbc))';
%             end 
%         end
        
%         function w = set_w(obj, K, w0)
%             if nargin < 3
%                 w0 = randn(obj.value.apx_function.N, 1);
%             end
%             [X,Y] = meshgrid(-2:.1:2, -2:.1:2);
%             mesh_size = size(X, 1);
%             Z = zeros(mesh_size, mesh_size);
%             for itr1 = 1 : mesh_size
%                for itr2 = 1 : mesh_size
%                     Z(itr1, itr2) = obj.cumulative_reward(K, [X(itr1,itr2);Y(itr1,itr2)]);
%                end
%             end
%             state = [vec(X), vec(Y)];
% %             Z = -1*ones(size(Z));
%             target = vec(Z);
% 
%             options = optimoptions('lsqnonlin','Display','iter','SpecifyObjectiveGradient',true);
%             w = lsqnonlin(@(w)obj.apx_cost_function(state, target, w),w0,[],[],options);
% % %             
% % %             [X,Y] = meshgrid(-0.5:.1:0.5, -2:.4:2);
% % %             mesh_size = size(X, 1);
% % %             Z = zeros(mesh_size, mesh_size);
% % %             for itr1 = 1 : mesh_size
% % %                for itr2 = 1 :mesh_size
% % %                     Z(itr1, itr2) = obj.value.est_value([X(itr1,itr2),Y(itr1,itr2)], w);
% % %                end
% % %             end
% % %             mesh(X,Y,Z)
%             
%         end
    end
end

