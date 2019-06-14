classdef RL_state_feedback_and_observer_train < RL_state_feedback_train
    %UNTITLED3 このクラスの概要をここに記述
    %   詳細説明をここに記述    
    methods
        function [true_x_all, apx_x_all, u_mpc_all, u_rl_all, omega] = actor_critic(obj, ini)
            apx_nx = obj.model.apx_nx;
            true_nx = obj.model.true_nx;
            nu = obj.model.nu;
            ny = obj.model.ny;
            y_all = zeros(obj.sim_N, ny);
            apx_x_all = zeros(obj.sim_N, apx_nx);
            true_x_all = zeros(obj.sim_N, true_nx);
            u_mpc_all = zeros(obj.sim_N, nu);
            u_rl_all = zeros(obj.sim_N, nu);
            true_x_all(1, :) = ini';
            try
                y_all(1, :) = obj.model.true_sys.C*ini + obj.model.true_sys.D*0;
            catch
                % modelがtest_apx_model出ないとtrue_sysが存在しない
                y_all(1, :) = ini';
            end
            if apx_nx == true_nx
                apx_x_all(1, :) = ini';
            else
                apx_x_all(1, :) = zeros(1, obj.model.apx_nx);
            end
            % only first episode initialize
            theta = randn(obj.basis_N, 1);
            omega = zeros(obj.basis_N, nu);
            % model based control Gain
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            for episode = 1 : obj.max_episode
                % initialize every episode
                z_theta = zeros(obj.basis_N, 1);
                z_omega = zeros(obj.basis_N, nu);
                zeta = 1;
                r_all = 0;
                for itr = 1 : obj.sim_N-1
                    if itr ~= 1
                        apx_x_all(itr, :) = obj.model.observer(apx_x_all(itr-1, :), u_mpc_all(itr-1, :)+u_rl_all(itr-1, :), y_all(itr, :));
                        r = obj.reward(apx_x_all(itr, :), u_mpc_all(itr-1, :)+u_rl_all(itr-1, :)); % REWARD is got at previous time.
                        r_all = r_all+obj.gamma^(itr-1)*r;
                        V_k0 = obj.state_basis_func(apx_x_all(itr, :))'*theta;
                        V_k1 = obj.state_basis_func(apx_x_all(itr-1, :))'*theta;
                        delta = r + obj.gamma*V_k0 - V_k1;
                        z_theta = obj.gamma*obj.lambda_theta*z_theta + zeta*obj.state_basis_func(apx_x_all(itr-1, :));
                        e_k1 = (u_rl_all(itr-1, :)-obj.state_basis_func(apx_x_all(itr-1, :))'*omega)/obj.sigma2_pi*obj.state_basis_func(apx_x_all(itr-1, :));
%                             e_k1 = obj.state_basis_func(x_all(itr-1, :));
                        z_omega = obj.gamma*obj.lambda_omega*z_omega + zeta*e_k1;
                        theta = theta + obj.alpha*delta*z_theta;
                        omega = omega + obj.beta*delta*z_omega;
                        zeta = obj.gamma2*zeta;
%                             u_rl = obj.state_basis_func(x_all(itr, :))'*omega;
                    end
                    u_rl = obj.state_basis_func(apx_x_all(itr, :))'*omega + obj.sigma2_pi*randn(1, nu);
                    u_mbc = -K*apx_x_all(itr, :)';
%                     u_mbc = 0;
                    [ne_x, ne_y] = obj.model.dynamics(true_x_all(itr,:), u_mbc+u_rl);
                    true_x_all(itr+1, :) = ne_x';
                    y_all(itr+1, :) = ne_y';
                    u_mpc_all(itr, :) = u_mbc;
                    u_rl_all(itr, :) = u_rl;
                end
                subplot(3,1,1)
                plot(obj.t, apx_x_all)
                title(['Episode-',num2str(episode)])
                grid on
                subplot(3,1,2)
                cost_history(episode) = obj.cost(apx_x_all, u_mpc_all+u_rl_all);
                plot(nonzeros(cost_history),'-r')
                ylabel('Cost')
                subplot(3,1,3)
                reward_history(episode) = r_all;
                plot(nonzeros(reward_history),'-b')
                ylabel('Integral Reward')
                drawnow
            end
        end
    end
end

