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
        
        function ne_x = dynamics(obj, pre_x, pre_u)
            func = @(x, u) [
                    x(2); ...
                    (obj.g/obj.l*sin(x(1)) - obj.eta/(obj.M*obj.l^2)*x(2) + 1/(obj.M*obj.l^2)*u);
                    ];
            ne_x = obj.RK4(func, pre_x, pre_u);
%             ne_x = obj.Euler(func, pre_x, pre_u);
            y = ne_x;
            
%             [
%                 x(1) + Ts*x(2), ...
%                 x(2) + Ts*(obj.g/obj.l*sin(x(1)) - obj.eta/(obj.M*obj.l^2)*x(2) + 1/(obj.M*obj.l^2)*u);
%             ];
        end        

        function ne_x = apx_dynamics(obj, pre_x, pre_u)
            system = @(x, u) obj.A*x + obj.B*u;
            ne_x = obj.RK4(system,pre_x, pre_u);          
        end
        
        function u = control_law(obj, set)
            switch set.law
                case 'lqr'
                    u = -set.K*set.x;
            end
        end
    end

end

