classdef CRLMBC_test_model < RL_model
    % 摩擦有倒立振子のモデル
    
    properties
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
        l
        M
        g
        eta
    end
    
    methods
        function obj = CRLMBC_test_model(l, M, g, eta, Ts)
            obj.A = eye(2) + Ts*[0, 1;g/l, -eta/(M*l^2)];
            obj.B = Ts*[0; 1/(M*l^2)];
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
        
        function [ne_x, y] = dynamics(obj, pre_x, pre_u)
            ne_x = [
                    pre_x(1) + obj.Ts*pre_x(2);
                    pre_x(2) + obj.Ts*(obj.g/obj.l*sin(pre_x(1)) - obj.eta/(obj.M*obj.l^2)*pre_x(2) + 1/(obj.M*obj.l^2)*pre_u);
                    ];
            y = ne_x;
        end        

        function ne_x = approximate_dynamics(obj, pre_x, pre_u)
            ne_x = obj.A*pre_x + obj.B*pre_u;            
        end
        
    end
end

