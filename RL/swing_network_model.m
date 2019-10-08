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
        system_rect
        dynamics_all
        observe_all
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
            % rect_sys
            rect_sys = Retrofit.generate_rectifier(obj.sys_local, apx_environment);
            rect_sys = c2d(rect_sys, obj.Ts, 'foh');
%             rect_sys = rect_sys({'x'},:);
            rect_sys = rect_sys({'yhat','what'},:);
            obj.system_rect = [rect_sys.A, rect_sys.B; rect_sys.C, rect_sys.D];
            % all_sys
            sys_all = obj.sys_all({obj.port_y,obj.port_w,obj.port_v}, {obj.port_control,obj.port_d_L}); 
            sys_all = c2d(sys_all, obj.Ts, 'foh');
            obj.observe_all  = sys_all.C;
            obj.dynamics_all = [sys_all.A, sys_all.B];
        end
        
% %         function set_controlled_system(obj, Q, R)
% %             obj.net.controllers = {};
% %             obj.net.add_controller(obj.c_n, obj.apx_sys_env, Q, R);
% %             sys_controlled = obj.net.get_sys_controlled(obj.sys_all);
% %             sys_controlled_extract = sys_controlled({obj.port_y, obj.port_xhat, obj.port_w, obj.port_v}, {obj.port_d_L, obj.port_control});
% %             sys_controlled_extract = c2d(sys_controlled_extract, obj.Ts, 'foh');
% % %             x_dim = obj.local_nx + obj.env_nx + obj.rect_nx;
% % %             y_dim = obj.ny + obj.ny + obj.w + obj.v; % y, xhat, w, v
% % %             u_dim = 2 + obj.nu;% noise port, control input
% %             [a,b,c,d] = ssdata(sys_controlled_extract);
% %             obj.system = [a,b;c,d;];
% %         end
        
%         function [local_ne_x, env_ne_x, ne_ywv, rect_ne_x] = dynamics(obj, local_pre_x, env_pre_x, rect_pre_x, all_pre_ywv, u, d_L)
% %             all_pre_x = obj.L_and_E_to_all(local_pre_x, env_pre_x);
% %             all_ne_x  = obj.RK4(obj.system_ywv_con_u, all_pre_x, u);
% %             ne_ywv    = obj.observer_ywv_con_u(all_ne_x, u);
% %             rect_ne_x = obj.RK4(obj.system_rect, rect_pre_x, all_pre_ywv);
% %             [local_ne_x, env_ne_x] = obj.all_to_L_ans_E(all_ne_x);
%             all_pre_x = obj.L_and_E_to_all(local_pre_x, env_pre_x);
%             xy = obj.system*[all_pre_x;rect_pre_x; d_L ;u];
%             all_ne_x = xy(1 : obj.local_nx+obj.env_nx);
%             rect_ne_x = xy(obj.local_nx+obj.env_nx+1 : obj.local_nx+obj.env_nx+obj.rect_nx);
%             ne_ywv = xy(obj.local_nx+obj.env_nx+obj.rect_nx+1 : end);
%             [local_ne_x, env_ne_x] = obj.all_to_L_ans_E(all_ne_x);
%         end
%         
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
        
        function [ywv, ne_local_x, ne_env_x] = dynamics(obj, local_x, env_x, u, d_L)
            narginchk(3, inf)
            all_x = obj.L_and_E_to_all(local_x, env_x);
            if nargout == 1 
                ywv = obj.observe_all*all_x;
            else
                ywv = obj.observe_all*all_x; % a little void
                ne_all_x = obj.dynamics_all*[all_x; d_L ;u];
                [ne_local_x, ne_env_x] = obj.all_to_L_ans_E(ne_all_x);
            end
        end
        
        function [ne_rect_x, rect_yw] = rect_dynamics(obj, local_x, env_x, rect_x)
            all_x = obj.L_and_E_to_all(local_x, env_x);
            ywv = obj.observe_all*all_x;
            rect_xy = obj.system_rect*[rect_x;ywv];
            ne_rect_x = rect_xy(1:obj.rect_nx, :);
            rect_yw = rect_xy(obj.rect_nx+1:end, :);
        end
    end
end

