classdef CRLMBC_test_model < environment_model
    % 摩擦有倒立振子のモデル
    
    properties
        A
        B
        C
        D
        K
        l
        M
        g
        eta
        controlled_sys
    end
    
    methods
        function obj = CRLMBC_test_model(l, M, g, eta, Ts)
%             obj.A = eye(2) + Ts*[0, 1;g/l, -eta/(M*l^2)];
            obj.A = [0, 1;g/l, -eta/(M*l^2)];
%             obj.B = Ts*[0; 1/(M*l^2)];
            obj.B = [0; 1/(M*l^2)];
            obj.C = eye(2);
            obj.D = zeros(2, 1);
            obj.Ts = Ts;
            obj.true_nx = 2;
            obj.apx_nx = 2;
            obj.nu = 1;
            obj.ny = 2;
            obj.l = l;
            obj.M = M;
            obj.g = g;
            obj.eta = eta;
        end
        

        function set_controlled_sys(obj, K)
            func = @(x, u) [
                    x(2); ...
                    (obj.g/obj.l*sin(x(1)) - obj.eta/(obj.M*obj.l^2)*x(2) + 1/(obj.M*obj.l^2)*(u-K*x));
                    ];
            obj.controlled_sys = func;
        end
        
        function [ne_x, y] = control_dynamics(obj, pre_x, pre_u)
            ne_x = obj.RK4(obj.controlled_sys, pre_x, pre_u);
            y = ne_x;
        end    
        
        function [ne_x, y] = dynamics(obj, pre_x, pre_u)
            func = @(x, u) [
                    x(2); ...
                    (obj.g/obj.l*sin(x(1)) - obj.eta/(obj.M*obj.l^2)*x(2) + 1/(obj.M*obj.l^2)*u);
                    ];
            ne_x = func(pre_x, pre_u);
            y = ne_x; 
        end

        function [ne_x, y] = approximate_dynamics(obj, pre_x, pre_u)
            func = @(x, u) obj.A*x + obj.B*u;
            ne_x = obj.RK4(func, pre_x, pre_u);
            y = obj.C*pre_x + obj.D*pre_u;
        end
        
    end

end

