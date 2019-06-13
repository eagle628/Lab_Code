classdef RL_state_feedback_and_observer_train < RL_state_feedback_train
    %UNTITLED3 このクラスの概要をここに記述
    %   詳細説明をここに記述    
    methods
        function [y_all, u_mpc_all, u_rl_all, omega] = actor_critic(obj, ini)
            nx = size(obj.model.A, 1);
            nu = size(obj.model.B, 2);
            ny = size(obj.model.C, 1);
            y_all = zeros(obj.sim_N, ny);
            x_all = zeros(obj.sim_N, nx);
            u_mpc_all = zeros(obj.sim_N, nu);
            u_rl_all = zeros(obj.sim_N, nu);
            y_all(1, :) = ini';
            x_all(1, :) = ini';
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
                    if itr == 1
                        u_rl = 0;
                    else
%                         if sum(abs(theta) > 0.5) == 0
                            x_all(itr, :) = obj.model.observer(x_all(itr-1, :), u_mpc_all(itr-1, :)+u_rl_all(itr-1, :), y_all(itr, :));
                            r = obj.reward(x_all(itr, :), u_mpc_all(itr-1, :)+u_rl_all(itr-1, :)); % REWARD is got at previous time.
                            r_all = r_all+obj.gamma^(itr-1)*r;
                            V_k0 = obj.state_basis_func(x_all(itr, :))'*theta;
                            V_k1 = obj.state_basis_func(x_all(itr-1, :))'*theta;
                            delta = r + obj.gamma*V_k0 - V_k1;
                            z_theta = obj.gamma*obj.lambda_theta*z_theta + zeta*obj.state_basis_func(x_all(itr-1, :));
                            e_k1 = (u_rl_all(itr-1, :)-obj.state_basis_func(x_all(itr-1, :))'*omega)/obj.sigma2_pi*obj.state_basis_func(x_all(itr-1, :));
%                             e_k1 = obj.state_basis_func(x_all(itr-1, :));
                            z_omega = obj.gamma*obj.lambda_omega*z_omega + zeta*e_k1;
                            theta = theta + obj.alpha*delta*z_theta;
                            omega = omega + obj.beta*delta*z_omega;
                            zeta = obj.gamma2*zeta;
                            u_rl = obj.state_basis_func(x_all(itr, :))'*omega + obj.sigma2_pi*randn(nu, 1);
%                             u_rl = obj.state_basis_func(x_all(itr, :))'*omega;
%                         else
%                             r = -10;
%                             r_all = r_all+r;
%                             V_k0 = 0;
%                             V_k1 = obj.state_basis_func(x_all(itr-1, :))'*theta;
%                             delta = r + obj.gamma*V_k0 - V_k1;
%                             z_theta = obj.gamma*obj.lambda_theta*z_theta + zeta*obj.state_basis_func(x_all(itr-1, :));
%                             e_k1 = (u_rl_all(itr-1, :)-obj.state_basis_func(x_all(itr-1, :))'*omega)/obj.sigma2_pi*obj.state_basis_func(x_all(itr-1, :));
%                             z_omega = obj.gamma*obj.lambda_omega*z_omega + zeta*e_k1;
%                             theta = theta + obj.alpha*delta*z_theta;
%                             omega = omega + obj.beta*delta*z_omega;
% %                             zeta = obj.gamma*zeta;
% %                             u_rl = obj.state_basis_func(x_all(itr, :))'*omega + obj.sigma2_pi*randn(1);
%                             break;
%                         end
                    end
                    u_mbc = -K*x_all(itr, :)';
%                     u_mbc = 0;
                    ne_x = obj.model.dynamics(y_all(itr,:), u_mbc+u_rl);
                    y_all(itr+1, :) = ne_x';
                    u_mpc_all(itr, :) = u_mbc;
                    u_rl_all(itr, :) = u_rl;
                end
                subplot(3,1,1)
                plot(obj.t, y_all)
                title(['Episode-',num2str(episode)])
                grid on
                subplot(3,1,2)
                plot(episode, obj.cost(y_all, u_mpc_all+u_rl_all),'*r')
                ylabel('Cost')
                hold on
                subplot(3,1,3)
                plot(episode, r_all,'*b')
                ylabel('Integral Reward')
                hold on
                drawnow
            end
        end
    end
end

