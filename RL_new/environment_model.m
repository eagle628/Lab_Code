classdef environment_model < handle
    % Environment on Reinforcement learning
    % Args
    % Ts : time step
    % nx : state number
    % nu : input number
    % ny : output number
    % state0 : current state memory
    % state1 : neext state
    
    properties
        Ts
        nx
        nu
        ny
        state
    end
    
    methods(Abstract)
        dynamics(obj)
    end
    
    methods
        
%         function x_hat = observer(obj, pre_x_hat, u, y)
%             if isempty(obj.L)
%                 obj.set_observer_gain();
%             end
%             x_hat = obj.A*pre_x_hat' + obj.B*u' - obj.L*(y' - obj.C*obj.A*pre_x_hat' - obj.C*obj.B*u' - obj.D*u');
%         end
%         
%         function set_observer_gain(obj, max)
%             if nargin < 2
%                 max = 1;
%             end
%             Flag = 1;
%             while Flag
%                 test_L = randn(obj.apx_nx, obj.ny);
%                 stable = (eye(obj.apx_nx, obj.apx_nx) + test_L*obj.C)*obj.A;
%                 stable = eig(stable);
%                 Flag = sum(abs(stable) > max);
%             end
%             obj.L = test_L;
%         end

        function ne_x = RK4(obj, system, pre_x, pre_u)
            k1 = system(pre_x, pre_u);
            k2 = system(pre_x+obj.Ts/2*k1, pre_u);
            k3 = system(pre_x+obj.Ts/2*k2, pre_u);
            k4 = system(pre_x+obj.Ts*k3, pre_u);
            ne_x = pre_x+obj.Ts/6*(k1+2*k2+2*k3+k4);
        end
        
        function ne_x = Euler(obj, system, pre_x, pre_u)
            x_dot = system(pre_x, pre_u);
            ne_x = pre_x+obj.Ts*x_dot;
        end

        function initialize(obj, ini_state)
            obj.state = ini_state;
        end
    end
end

