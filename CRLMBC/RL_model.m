classdef RL_model < handle
    %UNTITLED4 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Abstract)
        A
        B
        C
        D
        L
        Ts
        true_nx
        apx_nx
        nu
        ny
    end
    
    methods(Abstract)
        dynamics(obj)
    end
    
    methods
        function [ne_x, y] = approximate_dynamics(obj, pre_x, pre_u)
            ne_x = obj.A*pre_x + obj.B*pre_u;
            y = obj.C*pre_x + obj.D:pre_u;
        end
        
        function x_hat = observer(obj, pre_x_hat, u, y)
            if isempty(obj.L)
                obj.set_observer_gain();
            end
            x_hat = obj.A*pre_x_hat' + obj.B*u' - obj.L*(y' - obj.C*obj.A*pre_x_hat' - obj.C*obj.B*u' - obj.D*u');
        end
        
        function set_observer_gain(obj, max)
            if nargin < 2
                max = 1;
            end
            Flag = 1;
            while Flag
                test_L = randn(obj.apx_nx, obj.ny);
                stable = (eye(obj.apx_nx, obj.apx_nx) + test_L*obj.C)*obj.A;
                stable = eig(stable);
                Flag = sum(abs(stable) > max);
            end
            obj.L = test_L;
        end

        function ne_x = RK4(obj, system, pre_x, pre_u)
            k1 = system(pre_x, pre_u);
            k2 = system(pre_x+obj.Ts/2*k1, pre_u);
            k3 = system(pre_x+obj.Ts/2*k2, pre_u);
            k4 = system(pre_x+obj.Ts*k2, pre_u);
            ne_x = pre_x+obj.Ts/6*(k1+2*k2+2*k3+k4);
        end
    end
end

