classdef Global_Target
    %UNTITLED2 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        controller
        obs_interpreter
        AE
        BE
        CE
        R
        S
        Ap
        Bp
        Cp
        Ts
    end
    
    methods
        function obj = Global_Target(model, controller, obs_interpreter)
            if nargin < 3
                 obs_interpreter = observation_accumulater(model.ny, 1);
            end
            env = model.RL_env_all_prime;
            env = obs_interpreter.connect(env);% add
            obj.AE = env.A;
            obj.BE = env(:, strcat('u_node',num2str(model.c_n))).B;
            obj.CE = env({'ywhat_set'},:).C;
            obj.R  = env(:, strcat('d_node',num2str(model.c_n))).B;
            obj.S  = env('y', :).C;
            obj.Ts = model.Ts;
            obj.controller = controller;
            local = model.sys_local_discrete;
%             [obj.Ap,obj.Bp,obj.Cp,~] = ssdata(local);
            obj.obs_interpreter = obs_interpreter;
            [obj.Ap,obj.Bp,obj.Cp] = obj.obs_interpreter.connect(local);
        end
        
        function [f, grad, y] = eval_func(obj, theta, ddd)
            obj.controller.set_params(theta);
            [Ak,Bk,Ck,Dk,dAk,dBk,dCk,dDk] = obj.controller.get_ss();

            Anew = [obj.AE+obj.BE*Dk*obj.CE, obj.BE*Ck; Bk*obj.CE, Ak];
            Bnew = [obj.R; tools.zeros(Ak, obj.R)];
            Cnew = [obj.S, tools.zeros(obj.S, Ak); Dk*obj.CE, Ck];
            sys = balreal(ss(Anew, Bnew, Cnew, [], obj.Ts));

            f = 0;
            grad = zeros(length(theta), 1);
            if length(ddd) == 3
                ddd = randn(ddd);
            end
            for iter1 = 1 : size(ddd, 3)
                [y,~,x] = lsim(sys, ddd(:,:,iter1), []);
                f = f + norm(y);
                if nargin == 2 
                    delta = sqrt(eps)^(1/3);
                    for iter2 = 1 : length(dAk)
                        one_hot = zeros(size(theta));
                        one_hot(itr1) = 1;
                        grad(iter2) = -(obj.eval_func(theta+delta*one_hot,ddd)-obj.eval_func(theta-delta*one_hot,ddd))/(2*delta);
                    end
                end
            end
            f = f/size(ddd, 3);
        end
        
        function [c, ceq, Gc, Gceq] = stable_con(obj, theta)
            obj.controller.set_params(theta);
            [Ak,Bk,Ck,Dk] = obj.controller.get_ss();
            A_all = [Ak,Bk*obj.Cp; obj.Bp*Ck, obj.Ap+obj.Bp*Dk*obj.Cp];
            Flag =  any(any(isnan(A_all))) || any(any(isinf(A_all)));
            if Flag
                pole = 100;
            else
                pole = eig(A_all);
            end
            c = max(abs(pole))- 1;
            ceq = [];
            if nargout > 3
                Gc = zeros(size(theta))';
                delta = sqrt(eps)^(1/3);
                for itr1 = 1 : length(theta)
                    one_hot = zeros(size(theta));
                    one_hot(itr1) = 1;
                    Gc(itr1) = -(obj.stable_con(theta+delta*one_hot)-obj.stable_con(theta-delta*one_hot))/(2*delta);
                end
                Gceq = [];
            end
        end
        
        function area_gard = finite_drivative(theta)
            step_sizes = [1e-8];
            opt_eval = eval_func(theta, ddd);
            area_grad = zeros(1, obj.controller.N, length(step_sizes));
            for ivar = 1 : obj.controller.N
                one_hot = zeros(1, obj.controller.N);
                one_hot(ivar) = 1;
                for k = 1 : length(step_sizes)
                    area_grad(1, ivar, k) = 1/step_sizes(k)*(eval_func(theta+one_hot*step_sizes(k), ddd)-opt_eval);
                end
            end 
        end
    end
    
    methods(Static)
        function out = PSOcon(f, g)
            if g >= 0
                out = g;
            else
                out = atan(f) - pi/2;
                out = 1e6*out;
            end
        end
    end
end

