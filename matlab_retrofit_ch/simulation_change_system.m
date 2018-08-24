function [y,v,w] = simulation_change_system(sys_ori,sys_add,u,t,ob_y,ID_in,add_initial,ob_xhat,con_x,ch_t)
    % If A single Controller
    % 途中でシステムが変更された場合
    if nargin < 10
        ch_t = round(t(end)/2);
    end
    if nargin < 9
        con_x = [];
    end
    if nargin < 8
        ob_xhat = [];
    end
    if nargin < 7
        add_initial = [0,0];
    end
    Ts = t(2)-t(1);
    ch_t = round(ch_t/Ts)-1;
    y_before = lsim(sys_ori(ob_y,ID_in),u(1:ch_t,:),t(1:ch_t));
    x = lsim(sys_ori({'x'},ID_in),u(1:ch_t,:),t(1:ch_t));
    x = [x(end,:),add_initial];% If Network add a single node 
    % Extend and Simple
    if ~isempty(ob_xhat)
        x_o_p = lsim(sys_ori(ob_xhat,ID_in),u(1:ch_t,:),t(1:ch_t));
        x = [x,x_o_p(end,:)];
    end
    % only Extend
    if ~isempty(con_x)
        x_c_s = lsim(sys_ori(con_x,ID_in),u(1:ch_t,:),t(1:ch_t));
        x = [x,x_c_s(end,:)];
    end
    u = u(ch_t:end,:);
    t = 0:Ts:Ts*(length(u)-1);
    y_after  = lsim(sys_add(ob_y,ID_in),u,t,x);
    if nargout > 1
        v  = lsim(sys_add('v_node1',ID_in),u,t,x);
        w  = lsim(sys_add('w_node1',ID_in),u,t,x);
    end
    y = [y_before;y_after(2:end,:)];
end