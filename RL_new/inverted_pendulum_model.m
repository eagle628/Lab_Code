classdef inverted_pendulum_model < environment_model
    % 摩擦有倒立振子のモデル
    
    properties
        A
        B
        l
        M
        g
        eta
        Q
        R
    end
    
    methods
        function obj = inverted_pendulum_model(l, M, g, eta, Ts)
            obj.A = eye(2) + Ts*[0, 1;g/l, -eta/(M*l^2)];
            obj.B = Ts*[0; 1/(M*l^2)];
            obj.Ts = Ts;
            obj.nx = 2;
            obj.nu = 1;
            obj.ny = 2;
            obj.l = l;
            obj.M = M;
            obj.g = g;
            obj.eta = eta;
            obj.Q = eye(2);
            obj.R = eye(1);
            obj.state = zeros(obj.nx, 1);
        end
        
        function [y, reward] = dynamics(obj, u)
            obj.state = [...
                obj.state(1) + obj.Ts*obj.state(2);
                obj.state(2) + obj.Ts*(obj.g/obj.l*sin(obj.state(1)) - obj.eta/(obj.M*obj.l^2)*obj.state(2) + 1/(obj.M*obj.l^2)*u);
                ];
            y = obj.observe(); 
            reward = -1/10*(y'*obj.Q*y + u'*obj.R*u);
        end
        
        function y = observe(obj)
           y = obj.state; 
        end
    end
end

