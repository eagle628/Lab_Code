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
        
        function [x_all, y_all, u_all, t, reward] = sim_lqrcontroller(obj, ini, Te, Q, R)
            if nargin < 5 || isempty(R)
                R = eye(obj.model.nu);
            end
            if nargin < 4 || isempty(Q)
                Q = eye(obj.model.local_nx);
            end
            t = (0:obj.model.Ts:Te)';
            obj.model.net.add_controller(obj.model.c_n, Q, R);
            sys_all = obj.model.net.get_sys_controlled(obj.model.sys_all);
            sys_all = sys_all([obj.model.port_y,obj.model.port_xhat], obj.model.port_d_L);
            u_len = size(sys_all.B, 2);
            ddd = zeros(length(t), u_len);
            sys_all = c2d(sys_all, obj.model.Ts, obj.model.discrete_type);      
            [xy_all, ~, x_all] = lsim(sys_all, ddd, t, ini);
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
                [~, ~, C2, ~] = ssdata(sys_controller('u', {'y', 'w', 'v'}));
                u_all(itr, :) = C2*xhat_all;
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
        
        function [x_all, y_all, t] = sim_original(obj, ini, Te)
            % ini : sys_all length cut off
            t = (0:obj.model.Ts:Te)';
            sys_all = obj.model.sys_all(obj.model.port_y, obj.model.port_d_L);
            sys_all = c2d(sys_all, obj.model.Ts, obj.model.discrete_type);
            u_len = size(sys_all.B, 2);
            ddd = zeros(length(t), u_len);
            [y_all, ~, x_all] = lsim(sys_all, ddd, t, ini(1:order(sys_all)));
            y_all = y_all';
            x_all = x_all';
        end
        
        function [x_all, y_all, u_all, t, reward] = sim_extendlqrcontroller(obj, ini, Te, apx_env_dim, Q, R)
            if nargin < 6 || isempty(R)
                R = eye(obj.model.nu);
            end
            if nargin < 5 || isempty(Q)
                Q = eye(obj.model.local_nx);
            end
            t = (0:obj.model.Ts:Te)';
            obj.model.net.add_controller(obj.model.c_n, balred(obj.model.sys_env,apx_env_dim), Q, R);
            sys_all = obj.model.net.get_sys_controlled(obj.model.sys_all);
            sys_all = sys_all([obj.model.port_y,obj.model.port_xhat,{'controller_x1'}], obj.model.port_d_L);
            sys_all = c2d(sys_all, obj.model.Ts, obj.model.discrete_type);
            
            u_len = size(sys_all.B, 2);
            ddd = zeros(length(t), u_len);
            pre_ini = zeros(order(sys_all), 1);
            pre_ini(1:length(ini)) = ini;
            [xy_all, ~, x_all] = lsim(sys_all, ddd, t, pre_ini);
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
                [~, ~, C2, ~] = ssdata(sys_controller('u', {'y', 'w', 'v'}));
                u_all(itr, :) = C2*xhat_all;
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
        
        function render(obj, t, x_all, y_all, reward_history, episode)
            subplot(3,1,1)
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
            subplot(3,1,2)
            plot(t, y_all(1:obj.model.nz, :));
            ylabel('$\hat{y}$','Interpreter','latex')
            grid on
            subplot(3,1,3)
            tmp = nonzeros(reward_history);
            plot(tmp,'-b')
            ylim([2*median(tmp), 0])
            ylabel('Culumative Reward')
            grid on
            drawnow
        end
    end
end

