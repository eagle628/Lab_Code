classdef CRLMBC_test_model < handle
    % CRMBCの数値検証用プログラム
    
    properties
        l
        M
        g
        eta
        Ts
        A
        B
        C
        D
        L
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
            obj.C = eye(2);
            obj.D = zeros(2, 1);
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
            nx = size(obj.A, 1);
            ny = size(obj.C, 1);
            Flag = 1;
            while Flag
                test_L = randn(nx, ny);
                stable = (eye(nx, nx) + test_L*obj.C)*obj.A;
                stable = eig(stable);
                Flag = sum(abs(stable) > max);
            end
            obj.L = test_L;
        end
    end
end

