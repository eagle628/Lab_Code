classdef RL_state_feedback_train < handle
%     CRLMBC_test
    %UNTITLED3 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Constant)
        Q = diag([10,1])
        R = 1
        lambda_theta = 0.99
        lambda_omega = 0.99
        alpha = 0.05
        beta = 0.001
        gamma = 0.9
        gamma2 = 0.999
        max_episode = 1e3
    end
    
    properties
        model
        sim_N
        t
        basis_N
        c
        sigma2_basis
        sigma2_pi
    end
    
    methods
        function obj = RL_state_feedback_train(model, Te, basis_N)
            obj.model = model;
            obj.sim_N = Te/model.Ts + 1;
            obj.t = (0:model.Ts:Te)';
            obj.basis_N = basis_N;
%             obj.c = randn(obj.basis_N, size(obj.model.A, 1));
            m = 0.2*(-5:5)';
            n = 0.8*(-5:5)';
            obj.c = [kron(m,ones(11,1)),repmat(n,11,1)]; 
            obj.sigma2_basis = ones(obj.basis_N, 1);
            obj.sigma2_pi = 0.01;
        end
        
        function J = cost(obj, x, u)
            J = 0;
            for itr = 1 : obj.sim_N
                J = J + x(itr, :)*obj.Q*x(itr, :)' + u(itr, :)*obj.R*u(itr, :)';
            end
        end
        
        function r = reward(obj, x, u)
            r = -1/10*(x*obj.Q*x' + u*obj.R*u');
        end
        
        function phi = state_basis_func(obj, x)
            % x is row vector
            phi = exp(-sum((x-obj.c).^2, 2)./(2*obj.sigma2_basis));
        end
        
        function [x_all] = sim_noncontroller(obj, ini)
            x_all = zeros(obj.sim_N, 2);
            x_all(1, :) = ini';
            for itr = 1 : obj.sim_N-1
                x_all(itr+1, :) = obj.model.dynamics(x_all(itr,:), 0);
            end
        end
        
        function [y_all] = sim_lqrcontroller(obj, ini)
            y_all = zeros(obj.sim_N, 2);
            y_all(1, :) = ini';
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            for itr = 1 : obj.sim_N-1
                y_all(itr+1, :) = obj.model.dynamics(y_all(itr,:), -K*y_all(itr,:)');
            end
        end
        
        function [x_all, u_mpc_all, u_rl_all, omega] = actor_critic(obj, ini)
            nx = size(obj.model.A, 1);
            nu = size(obj.model.B, 2);
            x_all = zeros(obj.sim_N, nx);
            u_mpc_all = zeros(obj.sim_N, nu);
            u_rl_all = zeros(obj.sim_N, nu);
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
                            omega = omega + obj.beta*delta.*z_omega;
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
                    ne_x = obj.model.dynamics(x_all(itr,:), u_mbc+u_rl);
                    x_all(itr+1, :) = ne_x';
                    u_mpc_all(itr, :) = u_mbc;
                    u_rl_all(itr, :) = u_rl;
                end
                subplot(3,1,1)
                plot(obj.t, x_all)
                title(['Episode-',num2str(episode)])
                grid on
                subplot(3,1,2)
                plot(episode, obj.cost(x_all, u_mpc_all+u_rl_all),'*r')
                ylabel('Cost')
                hold on
                subplot(3,1,3)
                plot(episode, r_all,'*b')
                ylabel('Integral Reward')
                hold on
                drawnow
            end
        end
        
        function x_all = sim(obj, ini, omega)
            x_all = zeros(obj.sim_N, 2);
            x_all(1, :) = ini';
            K = dlqr(obj.model.A, obj.model.B, obj.Q, obj.R);
            for itr = 1 : obj.sim_N-1
                if itr == 1
                    u_rl = 0;
                else
                    u_rl = obj.state_basis_func(x_all(itr, :))'*omega;
                end
                u_mbc = -K*x_all(itr, :)';
                ne_x = obj.model.dynamics(x_all(itr,:), u_mbc+u_rl);
                x_all(itr+1, :) = ne_x';
            end
        end
    end
end

