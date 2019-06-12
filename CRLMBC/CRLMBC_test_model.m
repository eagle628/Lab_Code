classdef CRLMBC_test_model
    % CRMBCの数値検証用プログラム
    
    properties
        l
        M
        g
        eta
        Ts
        A
        B
        sim_N
    end
    
    methods
        function obj = CRLMBC_test_model(l, M, g, eta, Ts)
            obj.l = l;
            obj.M = M;
            obj.g = g;
            obj.eta = eta;
            obj.Ts = Ts;
            obj.A = eye(2) + Ts*[0, 1;g/l, -eta/(M*l^2)];
            obj.B = Ts*[0; 1/(M*l^2)];
        end
        
        function ne_x = dynamics(obj, pre_x, pre_u)
            ne_x = [
                    pre_x(1) + obj.Ts*pre_x(2);
                    pre_x(2) + obj.Ts*(obj.g/obj.l*sin(pre_x(1)) - obj.eta/(obj.M*obj.l^2)*pre_x(2) + 1/(obj.M*obj.l^2)*pre_u);
                    ];
        end
        
        function ne_x = approximate_dynamics(obj, pre_x, pre_u)
            ne_x = obj.A*pre_x + obj.B*pre_u;
        end
        
    end
end

