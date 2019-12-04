classdef swing_network_model < environment_model
    % network reinforcement learning Environment
    % Args
    % net : network class (Ex. network_swing_simple)
    % c_n : interested subsystem node number
    % sys_local : interested subsystem SS
    % sys_local_discrete : c2d(sys_local)
    % sys_env : network system except c_n node system SS
    % apx_sys_env : approximate network system except c_n node system SS
    % sys_all : All network sysmet SS
    % rect_nx : Rectifier state number
    % local_nx : local system(sys_local) state number
    % env_nx : Environment system(sys_env) state number
    % port_v : local system input interconnection v name
    % port_w : local system output interconnection w name
    % port_y : local system output observation y name
    % port_xhat : rectifier state name
    % port_d_L : local system input diturbance
    % port_control : local system input control
    % RL_env_all : Environment for Rl Angent SS (from {port_control} to
    % {port_y,port_w,port_v}))
    % dynamics_A : RL_env_all.A
    % dynamics_B : RL_env_all.B
    % dynamics_C : RL_env_all.C {port_y,port_w}
    % evaluate_C : evaluation point matrix
    % Q : reward weight (state)
    % R : reward weight (input)
    % discrete_type : system dicretize type('zoh', or 'foh')
    % nz : evaluation size

    properties(Constant)
    end

    properties
        net
        c_n
        sys_local
        sys_local_discrete
        sys_env
        apx_sys_env
        sys_all
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
        discrete_type
        nz
    end

    methods
        function obj = swing_network_model(net, c_n, Ts, apx_environment, discrete_type)
            if nargin < 5
                discrete_type = 'zoh';
            end
            if nargin < 4 || isempty(apx_environment)
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
            obj.discrete_type = discrete_type;
            obj.sys_local_discrete = make_new_local(c2d(obj.sys_local,Ts,obj.discrete_type), apx_environment);
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
            RL_env_all = c2d(RL_env_all, Ts, obj.discrete_type);
            obj.RL_env_all = RL_env_all;
            [obj.dynamics_A,obj.dynamics_B,obj.dynamics_C,~] = ssdata(RL_env_all({'yhat','what'},:));
            [~,~,obj.evaluate_C,~] = ssdata(RL_env_all({'y'},:));
            obj.state = zeros(order(RL_env_all), 1);
            obj.nx = order(RL_env_all);
            obj.ny = size(obj.sys_local.OutputGroup.y,2)+size(obj.sys_local.OutputGroup.w, 2);
            obj.nu = size(obj.sys_local.InputGroup.u, 2);
            obj.nz = size(obj.sys_local.OutputGroup.y,2);
            % weight
            obj.Q = eye(obj.nz);
            obj.R = eye(obj.nu);
        end

        function [ywv, r] = dynamics(obj, u)
            obj.state = obj.dynamics_A*obj.state + obj.dynamics_B*u;
            ywv = observe(obj);
            z = evaluate(obj);
            r = obj.reward(z, u);
        end

        function r = reward(obj, z, u)
            r = -(z'*obj.Q*z + u'*obj.R*u);
        end

        function ywv = observe(obj)
            ywv = obj.dynamics_C*obj.state;
        end

        function z = evaluate(obj, data)
            narginchk(1, inf)
            if nargin == 2
                z = obj.evaluate_C*data;
                return
            end
            z = obj.evaluate_C*obj.state;
        end

        function update = constraint(obj, target, new_params, data)
            [a2,b2,c2,d2] = target.apx_function.get_ss(new_params);
%             AB = data.belief_sys;
%             a1 = AB(:, 1:size(AB, 1));
%             b1 = AB(:, size(AB, 1)+1:end);
%             Ak = [a1, tools.zeros(a1, a2); b2*a1, a2];
%             Bk = [b1; b2*b1];
%             Ck = [d2*a1 c2];
%             Dk = d2*b1;
%             [Ap,Bp,Cp,~]  = ssdata(obj.sys_local_discrete);
%             A_all = [Ak,Bk*Cp; Bp*Ck, Ap+Bp*Dk*Cp];
            [AL,BL,CL] = data.belief_sys.connect(obj.sys_local_discrete);
            A_all = [a2,b2*CL; BL*c2, AL+BL*d2*CL];
            Flag =  any(any(isnan(A_all))) || any(any(isinf(A_all)));
            if Flag
                update = ~Flag;
                return
            end
            pole = eig(A_all);
            update = max(abs(pole))< 1;
        end
    end
end

%% local
function sys = make_new_local(local, apx_env)
    if isempty(apx_env.A)
       sys = local({'y','w'},{'u'});
       return;
    end
    [AE,BE,CE,DE] = ssdata(apx_env);
    AL = local.A;
    Bu = local(:,'u').B;
    Bv = local(:,'v').B;
    Cy = local('y',:).C;
    Cw = local('w',:).C;
    Anew = [AE, BE*Cw; Bv*CE, Bv*DE*Cw+AL];
    Bnew = [tools.zeros(AE,Bu);Bu];
    Cnew = [tools.zeros(Cy,AE), Cy; tools.zeros(Cw,AE), Cw];
    sys = ss(Anew, Bnew, Cnew, [], local.Ts);
    sys.InputGroup.u = local.InputGroup.u;
    sys.OutputGroup = local({'y','w'},:).OutputGroup;
    sys = balreal(sys);
end
