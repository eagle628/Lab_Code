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
        
        function [x_all, y_all, t] = sim_lqrcontroller(obj, ini, Te, varargin)
            t = (0:obj.model.Ts:Te)';
            obj.model.net.add_controller(obj.model.c_n,obj.model.Q,obj.model.R);
            sys_all = obj.model.net.get_sys_controlled(obj.model.sys_all);
            sys_all = sys_all(obj.model.port_y, obj.model.port_d_L);
            sys_all = c2d(sys_all, obj.model.Ts, obj.model.discrete_type);
            obj.model.net.controllers = {};
            [y_all, ~, x_all] = lsim(sys_all, ddd, t, ini);
        end
        
        function render(obj, t, x_all, y_all, reward_history, episode, update_chance)
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
            subplot(3,1,3)
            plot(nonzeros(reward_history),'-b')
            ylabel('Culumative Reward')
            drawnow
            disp(strcat('Episode-',num2str(episode),' : value  constraint update times : ', num2str(obj.opt_value.counter) ,'/',num2str(update_chance)))
            disp(strcat('Episode-',num2str(episode),' : policy constraint update times : ', num2str(obj.opt_policy.counter) ,'/',num2str(update_chance)))
            timer = toc;
            fprintf('This epoch %f[s], Estimated time to finish:%f [h].\n',timer, timer*(obj.max_episode-episode)/3600)
        end
    end
end

