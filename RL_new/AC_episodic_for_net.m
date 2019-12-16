classdef  AC_episodic_for_net < AC_episodic
    % Episodic Actor Critic trainer
    
    properties
    end
    
    methods
        function obj = AC_episodic_for_net(model, opt_policy, opt_value, belief_sys)
            if nargin < 4
                belief_sys = [];
            end
            obj@AC_episodic(model, opt_policy, opt_value, belief_sys);
        end
        
        % overrid for noise
        function [x_all, y_all, rl_u_all, t, reward] = sim(obj, ini, Te, d)
            sim_N = Te/obj.model.Ts + 1;
            t = (0:obj.model.Ts:Te)';
            if nargin < 4 || isempty(d)
                d = zeros(obj.model.nd, sim_N);
            else
                if length(d) == 1
                    rng(d)
                    d = randn(obj.model.nd, length(t));
                end
            end
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
                [y_all(:, k+1), r] = obj.model.dynamics_prime(rl_u_all(:, k), d(:, k));
                x_all(:, k+1) = obj.model.state;
                obj.belief_sys.update_internal_state(y_all(:, k), rl_u_all(:, k));
                belief_state =  obj.belief_sys.estimate(y_all(:, k+1));
                reward = reward + obj.gamma^(k-1)*r;
            end
        end
        
        
        function [x_all, y_all, u_all, t, reward] = sim_lqrcontroller(obj, ini, Te, d, Q, R)
            if nargin < 6 || isempty(R)
                R = eye(obj.model.nu);
            end
            if nargin < 5 || isempty(Q)
                Q = eye(obj.model.local_nx);
            end
            t = (0:obj.model.Ts:Te)';
            if nargin < 4 || isempty(d)
                d = zeros(obj.model.nd, length(t));
            else
                if length(d) == 1
                    rng(d)
                    d = randn(obj.model.nd, length(t));
                end
            end
            obj.model.net.add_controller(obj.model.c_n, Q, R);
            sys_all = obj.model.net.get_sys_controlled(obj.model.sys_all);
            sys_all = sys_all([obj.model.port_y,obj.model.port_xhat], obj.model.port_d_L);
            sys_all = c2d(sys_all, obj.model.Ts, obj.model.discrete_type);      
            [xy_all, ~, x_all] = lsim(sys_all, d, t, ini);
            idx = 0;
            for port_name = obj.model.port_y
                idx = idx + length(sys_all.OutputGroup.(port_name{:}));
            end
            y_all = xy_all(:, 1:idx)';
            yidx = idx;
            u_all = zeros(obj.model.nu, length(t));
            for itr = 1 : length(obj.model.net.controllers)
                port_len = length(sys_all.OutputGroup.(obj.model.port_xhat{itr}));
                xhat_all = xy_all(:, idx+(1:port_len))';
                sys_all = c2d(obj.model.sys_all(:, obj.model.port_control{itr}), obj.model.Ts, obj.model.discrete_type);
                sys_controller = obj.model.net.controllers{itr}.sys;
                [~, ~, C2, D2] = ssdata(sys_controller('u', {'y', 'w', 'v'}));
                u_all(itr, :) = C2*xhat_all;
                % 
                lcoal_node_idx = obj.model.net.controllers{itr}.nodes;
                G_Lap = obj.model.net.get_Laplacian_out(lcoal_node_idx);
                Cv = -G_Lap(lcoal_node_idx, :)*obj.model.sys_all('theta',:).c;   % Flip interconnection affect 
                Cw = obj.model.sys_all('theta',:).c(lcoal_node_idx, :);
                IDX = [lcoal_node_idx(:)*2-1, lcoal_node_idx(:)*2]';
                Cy = obj.model.sys_all('x',:).c(IDX(:),:);
                C1 = [Cy; Cw; Cv];
                u_all(itr, :) = u_all(itr, :) + D2*C1*x_all(:, 1:obj.model.net.N*2)';
                %
                idx = idx + port_len;
            end
            reward = cellfun(@(y,u)obj.model.reward(y,u),...
                                    mat2cell(y_all, yidx, ones(1,size(y_all,2))),...
                                    mat2cell(u_all, obj.model.nu, ones(1,size(u_all,2))));
            tmp = 0:length(t)-1;
            tmp = arrayfun(@(x)obj.gamma^x,tmp);
            reward = dot(reward, tmp);
            x_all = x_all';
            obj.model.net.controllers = {};
        end
        
        function [x_all, y_all, t] = sim_original(obj, ini, Te, d)
            % ini : sys_all length cut off
            t = (0:obj.model.Ts:Te)';
            if nargin < 4 || isempty(d)
                d = zeros(obj.model.nd, length(t));
            else
                if length(d) == 1
                    rng(d)
                    d = randn(obj.model.nd, length(t));
                end
            end
            sys_all = obj.model.sys_all(obj.model.port_y, obj.model.port_d_L);
            sys_all = c2d(sys_all, obj.model.Ts, obj.model.discrete_type);
            [y_all, ~, x_all] = lsim(sys_all, d, t, ini(1:order(sys_all)));
            y_all = y_all';
            x_all = x_all';
        end
        
        function [x_all, y_all, u_all, t, reward] = sim_extendlqrcontroller(obj, ini, Te, apx_env_dim, d, Q, R)
            if nargin < 7 || isempty(R)
                R = eye(obj.model.nu);
            end
            if nargin < 6 || isempty(Q)
                Q = eye(obj.model.local_nx);
            end
            t = (0:obj.model.Ts:Te)';
            if nargin < 5 || isempty(d)
                d = zeros(obj.model.nd, length(t));
            else
                if size(d) == 1
                    rng(d)
                    d = randn(obj.model.nd, length(t));
                end
            end
            obj.model.net.add_controller(obj.model.c_n, balred(obj.model.sys_env,apx_env_dim), Q, R);
            sys_all = obj.model.net.get_sys_controlled(obj.model.sys_all);
            sys_all = sys_all([obj.model.port_y,obj.model.port_xhat,{'controller_x1'}], obj.model.port_d_L);
            sys_all = c2d(sys_all, obj.model.Ts, obj.model.discrete_type);
            
            pre_ini = zeros(order(sys_all), 1);
            pre_ini(1:length(ini)) = ini;
            [xy_all, ~, x_all] = lsim(sys_all, d, t, pre_ini);
%             y_all = y_all';
%             x_all = x_all';
%             u_all = [];
%             reward = [];
            idx = 0;
            for port_name = obj.model.port_y
                idx = idx + length(sys_all.OutputGroup.(port_name{:}));
            end
            y_all = xy_all(:, 1:idx)';
            yidx = idx;
            u_all = zeros(obj.model.nu, length(t));
            % Number of controller is Single.
            for itr = 1 : length(obj.model.net.controllers)
                port_len = length(sys_all.OutputGroup.(obj.model.port_xhat{itr})) + length(sys_all.OutputGroup.('controller_x1'));
                xhat_all = xy_all(:, idx+(1:port_len))';
                sys_all = c2d(obj.model.sys_all(:, obj.model.port_control{itr}), obj.model.Ts, obj.model.discrete_type);
                sys_controller = obj.model.net.controllers{itr}.sys;
                [~, ~, C2, D2] = ssdata(sys_controller('u', {'y', 'w', 'v'}));
                u_all(itr, :) = C2*xhat_all;
                % 
                lcoal_node_idx = obj.model.net.controllers{itr}.nodes;
                G_Lap = obj.model.net.get_Laplacian_out(lcoal_node_idx);
                Cv = -G_Lap(lcoal_node_idx, :)*obj.model.sys_all('theta',:).c;   % Flip interconnection affect 
                Cw = obj.model.sys_all('theta',:).c(lcoal_node_idx, :);
                IDX = [lcoal_node_idx(:)*2-1, lcoal_node_idx(:)*2]';
                Cy = obj.model.sys_all('x',:).c(IDX(:),:);
                C1 = [Cy; Cw; Cv];
                u_all(itr, :) = u_all(itr, :) + D2*C1*x_all(:, 1:obj.model.net.N*2)';
                %
                idx = idx + port_len;
            end
            reward = cellfun(@(y,u)obj.model.reward(y,u),...
                                    mat2cell(y_all, yidx, ones(1,size(y_all,2))),...
                                    mat2cell(u_all, obj.model.nu, ones(1,size(u_all,2))));
            tmp = 0:length(t)-1;
            tmp = arrayfun(@(x)obj.gamma^x,tmp);
            reward = dot(reward, tmp);
            x_all = x_all';
            obj.model.net.controllers = {};
        end
        
        function render(obj, t, x_all, y_all, reward_history, episode, delta)
            subplot(4,1,1)
            z = obj.model.evaluate(x_all);
            plot(t, z)
            title(['\fontsize{16}','Episode-',num2str(episode)])
            ylabel('y')
            grid on
            lim = max(max(abs(z),[],'omitnan'));
            if isempty(lim) || lim == 0
                lim = 1;
            end
            ylim([-lim, lim]);
            subplot(4,1,2)
            plot(t, y_all(1:obj.model.nz, :));
            ylabel('$\hat{y}$','Interpreter','latex')
            grid on
            subplot(4,1,3)
            tmp = nonzeros(reward_history);
            plot(tmp,'-b')
%             ylim([-1e4, 0])
            ylabel('Culumative Reward')
            grid on
            subplot(4,1,4)
            plot(episode, norm(delta),'*')
            hold on
            drawnow
        end
    end
end

