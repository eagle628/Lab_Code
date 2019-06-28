classdef swing_network_model < environment_model
    %UNTITLED6 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Constant)
    end
    
    properties
        net
        sys_local
        sys_env
        sys_all
        system_rect
        observer_rect
        system_ywv_d_L
        observer_ywv_d_L
        c_n
        A
        B
        rect_x
    end
    
    methods
        function obj = swing_network_model(net, c_n, Ts)
            obj.net = net;
            obj.c_n = c_n;
            sys_all = obj.net.get_sys();
            for nnn = 1 :  obj.net.N
                sys_all = net1.add_io(sys_all, nnn, strcat('node',num2str(nnn)));
            end
            obj.sys_all = sys_all;
            [obj.sys_local, obj.sys_env] = net.get_sys_local(obj.c_n);
            sys_design = obj.sys_local('y', 'u');
            obj.A = sys_design.A;
            obj.B = sys_design.B;
            obj.Ts = Ts;
            obj.true_nx = 2*obj.net.N-2;
            obj.apx_nx = 2*obj.net.N-2;
            obj.rect_x = 2;
            obj.nu = 1;
            obj.ny = 1;
            
            rect_sys = Retrofit.generate_rectifier(sys_local, ss(0));
            [a_rect,b_rect,c_rect,d_rect] = ssdata(rect_sys({'x'},{'y','w','v'}));
            obj.system_rect = @(x, u) a_rect*x + b_rect*u;
            obj.observer_rect = @(x, i) c_rect*x + d_rect*u;

            sys_continous_d_L_to_ywv = sys_all({ob_y_p, ob_w_p, ob_v_p}, {ID_in_p, control_in_p});

            [a_ywv_d_L, b_ywv_d_L, c_ywv_d_L, d_ywv_d_L] = ssdata(sys_continous_d_L_to_ywv);
            obj.system_ywv_d_L = @(x, u) a_ywv_d_L*x + b_ywv_d_L*u;
            obj.observer_yvw_d_L = @(x, u) c_ywv_d_L*x + d_ywv_d_L*u;       
        end
        
        function [all_ne_x, all_ne_ywv, rect_ne_x] = dynamics(obj, all_pre_x, rect_pre_x, all_pre_y, u)
            all_ne_x  = obj.RK4(obj.system_ywv_d_L, all_pre_x, u);
            all_ne_ywv  = obj.RK4(obj.observer_ywv_d_L, all_ne_x, u);
            rect_ne_x = obj.RK4(obj.system_rect, rect_pre_x, all_pre_y);
        end
        
        function u = control_law(obj, set)
            u = (set.K*set.rect_x_all(k, :)')' - (set.K*set.ywv_all(k,1:2)')';
        end
    end
end

