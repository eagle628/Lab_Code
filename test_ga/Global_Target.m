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
                    for iter2 = 1 : length(dAk)
                        dAnew = [obj.AE+obj.BE*dDk{iter2}*obj.CE, obj.BE*dCk{iter2}; dBk{iter2}*obj.CE, dAk{iter2}];
                        dBnew = [obj.R; tools.zeros(Ak, obj.R)];
                        dCnew = [obj.S, tools.zeros(obj.S, Ak); dDk{iter2}*obj.CE, dCk{iter2}];
                        dsys = balrelal(ss(Anew, [dAnew,dBnew], Cnew, [dCnew,tools.zeros(dCnew, dBnew)], obj.Ts));
                        dy = lsim(dsys, [x, ddd(:,:,iter1)], []);
                        grad(iter2) = grad(iter2) + sum(sum(y.*dy));
                    end
                end
            end
            f = f/size(ddd, 3);
        end
        
        function [c, ceq] = stable_con(obj, theta)
            obj.controller.set_params(theta);
            [Ak,Bk,Ck,Dk] = obj.controller.get_ss();
            A_all = [Ak,Bk*obj.Cp; obj.Bp*Ck, obj.Ap+obj.Bp*Dk*obj.Cp];
            pole = eig(A_all);
            c= max(abs(pole))- 1;
            ceq = [];
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

