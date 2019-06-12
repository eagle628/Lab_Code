classdef LQR_IRL_VI_alg_train < handle
    %UNTITLED5 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Constant)
       S = eye(2)
       R = eye(1)
       gamma = 0.01
       max_episode = 1e3
    end
    
    properties
        model
        t
        sim_N
    end
    
    methods
        function obj = LQR_IRL_VI_alg_train(model, Te)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
        end
        
        function x_all = train(obj, ini)
            P = zeros(2,2);
            x_all = zeros(obj.sim_N, 2);
            x_all(1, :) =  ini';
            u_all = zeros(obj.sim_N, 1);
            K = zeros(1, 2);
            for episode = 1 : obj.max_episode
                for itr1 = 1 : obj.sim_N-1
                    u_all(itr1, :) = K*x_all(itr1, :)';
                    x_all(itr1+1, :) = obj.model.dynamics(x_all(itr1, :), u_all(itr1, :));
                end
                [P, K] = obj.updator(x_all, u_all, P);
                plot(x_all)
                title(['Episode-',num2str(episode)])
                drawnow
            end
        end
        
        function [P1, K] = updator(obj, x, u, P0)
            P1 = obj.value_update(x, u, 0, obj.t(end), P0);
            K  = obj.policy_update(P1);
        end
        
        function next_P = value_update(obj, x, mu, t, T, pre_P)
            func = @(tau, x, mu) exp(-obj.gamma*(tau-t))*(x'*obj.S*x + mu'*obj.R*mu);
            xT_next_P_x = 0;
            series_t = t:obj.model.Ts:t+T-obj.model.Ts; 
            for itr1 = 1 : length(series_t)-1
                xT_next_P_x = xT_next_P_x + obj.model.Ts/2 * (func(series_t(itr1), x(itr1, :)', mu(itr1, :)') + func(series_t(itr1+1), x(itr1+1, :)', mu(itr1+1, :)'));
            end
            xT_next_P_x = xT_next_P_x + exp(-obj.gamma*T)*x(end, :)*pre_P*x(end, :)';
            next_P = (x(1, :)'*x(1, :))\(x(1, :)'*xT_next_P_x*x(1, :))/(x(1, :)'*x(1, :));
        end
        
        function K = policy_update(obj, P)
            K = -obj.R\obj.model.B'*P;
        end
    end
end


%% local
