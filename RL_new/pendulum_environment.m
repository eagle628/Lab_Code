classdef pendulum_environment < environment_model
    % 摩擦有振子のモデル
    % A : dynamics matrix(discrete)
    % B : Input matrix(discrete)
    % Q,R : reward weight
    
    properties
        A
        B
        K
        l
        M
        g
        eta
        Q
        R
    end
    
    methods
        function obj = pendulum_environment(l, M, g, eta, Ts)
            obj.A = eye(2) + Ts*[0, 1;-g/l, -eta/(M*l^2)];
            obj.B = Ts*[0; 1/(M*l^2)];
            obj.Ts = Ts;
            obj.nx = 2;
            obj.nu = 1;
            obj.ny = 2;
            obj.l = l;
            obj.M = M;
            obj.g = g;
            obj.eta = eta;
            obj.Q = diag([10,1]);
            obj.R = eye(1);
            obj.state = zeros(obj.nx, 1);
        end
        
        function [y, reward] = dynamics(obj, u)
            obj.state = [...
                obj.state(1) + obj.Ts*obj.state(2);
                obj.state(2) + obj.Ts*(-obj.g/obj.l*sin(obj.state(1)) - obj.eta/(obj.M*obj.l^2)*obj.state(2) + 1/(obj.M*obj.l^2)*u);
                ];
            y = observe(obj);
            reward = -1/10*(y'*obj.Q*y + u'*obj.R*u);
        end
        
        function y = observe(obj)
           y = obj.state; 
        end
        
    end

end

