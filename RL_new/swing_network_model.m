classdef swing_network_model < environment_model
    %UNTITLED6 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Constant)
    end
    
    properties
        net
        sys_local
        sys_local_discrete
        sys_env
        apx_sys_env
        sys_all
        c_n
        A
        B
        rect_nx
        local_nx
        env_nx
        port_v
        port_w
        port_y
        port_xhat
        port_d_L
        port_control
        RL_env_all
        dynamics_A
        dynamics_B
        dynamics_C
        evaluate_C
        Q
        R
    end
    
    methods
        function obj = swing_network_model(net, c_n, Ts, apx_environment)
            if nargin < 4
                apx_environment = ss(0);
            end
            obj.net = net;
            obj.c_n = c_n;
            obj.apx_sys_env = apx_environment;
            sys_all = obj.net.get_sys();
            for nnn = 1 :  obj.net.N
                sys_all = net.add_io(sys_all, nnn, strcat('node',num2str(nnn)));
            end
            obj.sys_all = sys_all;
            [obj.sys_local, obj.sys_env] = net.get_sys_local(obj.c_n);
            discrete_mode = 'zoh';
            obj.sys_local_discrete = c2d(obj.sys_local({'y','w'},{'u'}), Ts,discrete_mode);
            sys_design = Retrofit.generate_fb(obj.sys_local, apx_environment);
            sys_design = c2d(sys_design('y','u'),Ts,discrete_mode);
            obj.A = sys_design.A;
            obj.B = sys_design.B;
            obj.Ts = Ts;
            obj.rect_nx = length(c_n)*2 + order(apx_environment); % retifier dim
            obj.local_nx = order(obj.sys_local); % local dim
            obj.env_nx = order(obj.sys_env);
            % port name
            obj.port_v  = {strcat('v_node',num2str(c_n))};
            obj.port_w  = {strcat('w_node',num2str(c_n))};
            obj.port_y  = {strcat('y_node',num2str(c_n))};
            obj.port_xhat = {strcat('xhat_controlled1')};
            obj.port_d_L = {strcat('d_node',num2str(c_n))};
            obj.port_control = {strcat('u_node',num2str(c_n))};
            % rect_sys
            rect_sys = Retrofit.generate_rectifier(obj.sys_local, apx_environment);
            % all_sys
            sys_all = obj.sys_all({obj.port_y,obj.port_w,obj.port_v}, {obj.port_control}); 
            [A2,B2,C2,D2] = ssdata(rect_sys);
            [A1,B1,C1,D1] = ssdata(sys_all);
            Ak = [A1, tools.zeros(A1, A2); B2*C1, A2];
            Bk = [B1; B2*D1];
            Ck = [D2*C1 C2];
            Dk = D2*D1;
            RL_env_all = ss(Ak,Bk,Ck,Dk);
            RL_env_all.OutputGroup = rect_sys.OutputGroup;
            RL_env_all.InputGroup = sys_all.InputGroup;
            RL_env_all = c2d(RL_env_all, Ts, discrete_mode);
            obj.RL_env_all = RL_env_all;
            [obj.dynamics_A,obj.dynamics_B,obj.dynamics_C,~] = ssdata(RL_env_all({'yhat','what'},:));
            [~,~,obj.evaluate_C,~] = ssdata(RL_env_all({'y'},:));
            obj.state = zeros(order(RL_env_all), 1);
            obj.nx = order(RL_env_all); 
            obj.ny = size(obj.sys_local.OutputGroup.y,2)+size(obj.sys_local.OutputGroup.w, 2);
            obj.nu = size(obj.sys_local.InputGroup.u, 2);
            % weight
            obj.Q = 1;
            obj.R = 1;
        end
        
        function [ywv, reward] = dynamics(obj, u)
            obj.state = obj.dynamics_A*obj.state + obj.dynamics_B*u;
            ywv = observe(obj);
            z = evaluate(obj);
            reward = z(2)'*obj.Q*z(2) + u'*obj.R*u;
        end
        
        function ywv = observe(obj)
            ywv = obj.dynamics_C*obj.state;
        end
        
        function z = evaluate(obj, data)
            narginchk(1, inf)
            if nargin == 2
                data = mat2cell(data, size(data,1),ones(1, size(data,2)));
                z = cellfun(@(x)obj.evaluate_C*x, data, 'UniformOutput', false);
                z = cell2mat(z);
                return
            end
            z = obj.evaluate_C*obj.state; 
        end
    end
end

