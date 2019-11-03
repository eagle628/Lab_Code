classdef  AC_episodic_for_net < AC_episodic
    % Episodic Actor Critic trainer
    
    properties
    end
    
    methods
        function obj = AC_episodic_for_net(model, opt_policy, opt_value, belief_sys)
            if nargin < 4
                belief_sys = [zeros(model.ny), eye(model.ny)];
            end
            obj@AC_episodic(model, opt_policy, opt_value, belief_sys);
            obj.gamma = 0.99;
            obj.max_episode = 10000;
            obj.snapshot = 1000;
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
            sys_all = sys_all(obj.model.port_y, obj.model.port_d_L);
            u_len = size(sys_all.B, 2);
            ddd = zeros(length(t), u_len);
            sys_all = c2d(sys_all, obj.model.Ts, obj.model.discrete_type);
            obj.model.net.controllers = {};
            [y_all, ~, x_all] = lsim(sys_all, ddd, t, ini);
            u_all = [];
            reward = [];
        end
        
        function [x_all, y_all, t] = sim_original(obj, ini, Te)
            % ini : sys_all length cut off
            t = (0:obj.model.Ts:Te)';
            sys_all = obj.model.sys_all(obj.model.port_y, obj.model.port_d_L);
            sys_all = c2d(sys_all, obj.model.Ts, obj.model.discrete_type);
            u_len = size(sys_all.B, 2);
            ddd = zeros(length(t), u_len);
            [y_all, ~, x_all] = lsim(sys_all, ddd, t, ini(1:order(sys_all)));
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
            sys_all = sys_all(obj.model.port_y, obj.model.port_d_L);
            sys_all = c2d(sys_all, obj.model.Ts, obj.model.discrete_type);
            obj.model.net.controllers = {};
            u_len = size(sys_all.B, 2);
            ddd = zeros(length(t), u_len);
            pre_ini = zeros(order(sys_all), 1);
            pre_ini(1:length(ini)) = ini;
            [y_all, ~, x_all] = lsim(sys_all, ddd, t, pre_ini);
            u_all = [];
            reward = [];
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

