classdef swing_network_model < environment_model
    %UNTITLED6 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Constant)
    end
    
    properties
        net
        sys_local
        sys_env
        apx_sys_env
        sys_all
%         system_rect
%         observer_rect
%         system_ywv_con_u
%         observer_ywv_con_u
        system
        c_n
        A
        B
        rect_nx
        local_nx
        env_nx
        nw
        nv
        port_v
        port_w
        port_y
        port_xhat
        port_d_L
        port_control
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
            sys_design = Retrofit.generate_fb(obj.sys_local, apx_environment);
            sys_design = sys_design('y','u');
            obj.A = sys_design.A;
            obj.B = sys_design.B;
            obj.Ts = Ts;
            obj.apx_nx = order(apx_environment); % apx environmetn dim
            obj.nw = length(c_n);
            obj.nv = length(c_n);
            obj.rect_nx = length(c_n)*2 + obj.apx_nx; % retifier dim
            obj.local_nx = order(obj.sys_local); % local dim
            obj.env_nx = order(obj.sys_env);
            obj.true_nx = 2*obj.net.N-obj.local_nx; % true environment dim
            obj.ny = size(obj.sys_local.OutputGroup.y, 2);
            obj.nu = size(obj.sys_local.InputGroup.u, 2);
            
            % port name
            obj.port_v  = {strcat('v_node',num2str(c_n))};
            obj.port_w  = {strcat('w_node',num2str(c_n))};
            obj.port_y  = {strcat('y_node',num2str(c_n))};
            obj.port_xhat = {strcat('xhat_controlled1')};
            obj.port_d_L = {strcat('d_node',num2str(c_n))};
            obj.port_control = {strcat('u_node',num2str(c_n))};
            % extract system
% %             rect_sys = Retrofit.generate_rectifier(obj.sys_local, apx_environment);
% %             [a_rect,b_rect,c_rect,d_rect] = ssdata(rect_sys({'x'},{'y','w','v'}));
% %             obj.system_rect = @(x, u) a_rect*x + b_rect*u;
% %             obj.observer_rect = @(x, i) c_rect*x + d_rect*u;
% % 
% %             sys_continous_con_u_to_ywv = sys_all({obj.port_y, obj.port_w, obj.port_v}, {obj.port_control});
% % 
% %             [a_ywv_con_u, b_ywv_con_u, c_ywv_con_u, d_ywv_con_u] = ssdata(sys_continous_con_u_to_ywv);
% %             obj.system_ywv_con_u = @(x, u) a_ywv_con_u*x + b_ywv_con_u*u;
% %             obj.observer_ywv_con_u = @(x, u) c_ywv_con_u*x + d_ywv_con_u*u;       
        end
        
        function set_controlled_system(obj, Q, R)
            obj.net.add_controller(obj.c_n, obj.apx_sys_env, Q, R);
            sys_controlled = obj.net.get_sys_controlled(obj.sys_all);
            sys_controlled_extract = sys_controlled({obj.port_y, obj.port_xhat, obj.port_w, obj.port_v}, {obj.port_d_L, obj.port_control});
            sys_controlled_extract = c2d(sys_controlled_extract, obj.Ts, 'foh');
%             x_dim = obj.local_nx + obj.env_nx + obj.rect_nx;
%             y_dim = obj.ny + obj.ny + obj.w + obj.v; % y, xhat, w, v
%             u_dim = 2 + obj.nu;% noise port, control input
            [a,b,c,d] = ssdata(sys_controlled_extract);
            obj.system = [a,b;c,d;];
        end
        
        function [local_ne_x, env_ne_x, ne_ywv, rect_ne_x] = dynamics(obj, local_pre_x, env_pre_x, rect_pre_x, all_pre_ywv, u)
%             all_pre_x = obj.L_and_E_to_all(local_pre_x, env_pre_x);
%             all_ne_x  = obj.RK4(obj.system_ywv_con_u, all_pre_x, u);
%             ne_ywv    = obj.observer_ywv_con_u(all_ne_x, u);
%             rect_ne_x = obj.RK4(obj.system_rect, rect_pre_x, all_pre_ywv);
%             [local_ne_x, env_ne_x] = obj.all_to_L_ans_E(all_ne_x);
            all_pre_x = obj.L_and_E_to_all(local_pre_x, env_pre_x);
            xy = obj.system*[all_pre_x;rect_pre_x;0;0;u];
            all_ne_x = xy(1 : obj.local_nx+obj.env_nx);
            rect_ne_x = xy(obj.local_nx+obj.env_nx+1 : obj.local_nx+obj.env_nx+obj.rect_nx);
            ne_ywv = xy(obj.local_nx+obj.env_nx+obj.rect_nx+1 : end);
            [local_ne_x, env_ne_x] = obj.all_to_L_ans_E(all_ne_x);
        end
        
        function all_x = L_and_E_to_all(obj, local_x, env_x)
            lidx = obj.generate_local_idx();
            all_x = zeros(2*obj.net.N, 1);
            all_x(lidx) = local_x;
            all_x(~lidx) = env_x;
        end
        
        function [local_x, env_x] = all_to_L_ans_E(obj, all_x)
            lidx = obj.generate_local_idx();
            local_x = all_x(lidx);
            env_x = all_x(~lidx);
        end
        
        function lidx = generate_local_idx(obj)
            lidx = false(2*obj.net.N, 1);
            for nidx = obj.c_n
               lidx(2*nidx-1) = true;
               lidx(2*nidx)   = true;
            end
        end
    end
end

