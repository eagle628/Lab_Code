%% Function Note
% INPUT
% Case : When the system is changed in the middel time.
% sys_ori : original system
% sys_ch : changed system
% ob_y : Observation point (Cell type)
% noise_point : Noise Insert point (Cell type)
% condition : Add Node -> Node Initial (EX : Add 1 Node is [0,0])
%             Add Controller -> Output Rectifer State and Extend Controller initial
%             Remove Controller -> i'th Controller is Removed ( scalar number "Not 1")
% ch_t : Change System time ( Not Step ,Real Time is)

% OUTPUT
% y : mesurement output
% xhat : output rectifier states
% v : interconnection signal (G_E output)
% vhat_con : prediction v in controller (????)
% w : ibnterconnection signal (G_E input)

% ATTENSION
% Add Controller append on last.
%% Main Function
function [ y, xhat, v, vhat_con, w] = sim_change_system(number, sys_ori, sys_ch, ob_y_p, noise_point, u, t, condition, ch_t)
    if nargin < 9 || isempty(ch_t) 
        ch_t = round(t(end)/2);
    end
    if nargin < 8
        condition = [];
    end
    % Simulation Time
    Ts = t(2)-t(1);
    ch_t = round(ch_t/Ts)+1;
    % System Information
    controller1 = sys_ori.UserData.controllers;
    controller2 = sys_ch.UserData.controllers;
    Node1 = sys_ori.UserData.Nodes;
    Node2 = sys_ch.UserData.Nodes;
    controller_pattern1 = sys_ori.UserData.controllertype;
    controller_pattern2 = sys_ch.UserData.controllertype;
    control_nodes1 = sys_ori.UserData.control_nodes;
    control_nodes2 = sys_ch.UserData.control_nodes;
    % Add Node
    if (controller1 <= controller2) && (Node1 < Node2)
        y_before = lsim(sys_ori(ob_y_p,noise_point),u(1:ch_t,:),t(1:ch_t));
        x = lsim(sys_ori({'x'},noise_point),u(1:ch_t,:),t(1:ch_t));
        if isempty(condition)
            condition = [0,0];
            %condition = [0,x(end,2)];
        end
        condition = kron(ones(1,Node2-Node1),condition);
        x = [x(end,:),condition];
        % Controller State
        for itr = 1 : controller1
            % Extend and Simple
            x_o_p = lsim( sys_ori(strcat('xhat_controlled',num2str(itr)), noise_point), u(1:ch_t,:), t(1:ch_t));
            x = [x,x_o_p(end,:)];
            if itr == number && nargout > 1
                xhat_before = x_o_p;
            end
            % only Extend
            if controller_pattern1(itr) == 'e'
                x_c_s = lsim(sys_ori( strcat('controller_x',num2str(itr)), noise_point) ,u(1:ch_t,:), t(1:ch_t));
                x = [x,x_c_s(end,:)];
            end
        end
        if controller1 < controller2
            x = [x, 0,0];
        end
        u_after = u(ch_t:end,:);
        t_after = t(ch_t:end,:) - t(ch_t);
        y_after = lsim( sys_ch(ob_y_p,noise_point), u_after, t_after, x);
        y = [y_before;y_after(2:end,:)];
        if nargout > 1
            xhat_after = lsim( sys_ch(strcat('xhat_controlled',num2str(number)), noise_point), u_after, t_after, x);
            xhat = [xhat_before;xhat_after(2:end,:)];
        end
        if nargout > 2
            v_before = lsim( sys_ori(strcat('v_node',num2str(control_nodes1{number})),noise_point),u(1:ch_t,:),t(1:ch_t));
            v_after = lsim( sys_ch(strcat('v_node',num2str(control_nodes1{number})),noise_point), u_after, t_after, x);
            v = [v_before;v_after(2:end,:);];
        end
        if nargout > 3
            vhat_con_before = lsim( sys_ori(strcat('vhat_controlled',num2str(number)),noise_point),u(1:ch_t,:),t(1:ch_t));
            vhat_con_after = lsim( sys_ch(strcat('vhat_controlled',num2str(number)),noise_point), u_after, t_after, x);
            vhat_con = [vhat_con_before;vhat_con_after(2:end,:);];
        end
        if nargout > 4
            w_before = lsim( sys_ori(strcat('w_node',num2str(control_nodes1{number})),noise_point),u(1:ch_t,:),t(1:ch_t));
            w_after = lsim( sys_ch(strcat('w_node',num2str(control_nodes1{number})),noise_point), u_after, t_after, x);
            w = [w_before;w_after(2:end,:);];
        end
    end
    % Add Controller ( Single Add Controller )
    if (Node1 == Node2) && (controller1 < controller2)
        % Add Controller
        y_before = lsim(sys_ori(ob_y_p,noise_point),u(1:ch_t,:),t(1:ch_t));
        x = lsim(sys_ori({'x'},noise_point),u(1:ch_t,:),t(1:ch_t));
        % Controller State
        for itr = 1 : controller1
            % Extend and Simple
            x_o_p = lsim(sys_ori(strcat('xhat_controlled',num2str(itr)), noise_point), u(1:ch_t,:), t(1:ch_t));
            x = [x(end,:),x_o_p(end,:)];
            if itr == number && nargout > 1
                xhat_before = x_o_p;
            end
            % only Extend
            if controller_pattern1(itr) == 'e'
                x_c_s = lsim(sys_ori( strcat('controller_x',num2str(itr)), noise_point), u(1:ch_t,:), t(1:ch_t));
                x = [x,x_c_s(end,:)];
            end
        end
        % add controller Sate
        add_control_nodes = control_nodes2{end};
        % Simple And Extend
            % output rectifer's state can measure or calucaulate.
%             for itr = 1 : numel(add_control_nodes)
%                 x = [ x, x(2*add_control_nodes(itr)-1:2*add_control_nodes(itr))];
%             end
        % It is natual that Both RetroFit Controller initial is [ 0, 0]. 
        x = [ x, 0, 0];
        % Extend only
        if controller_pattern2(end) == 'e'
            dim = numel(getfield( sys_ch.OutputGroup, strcat('controller_x',num2str(controller2))));
            if isempty(condition)
                condition = 0;
                x = [x,ones(1,dim)*condition];
            else
                if numel(condition) == dim
                    x = [x,condition];
                else
                    x = [x,ones(1,dim)*condition(1)];
                end
            end
        end
        u_after = u(ch_t:end,:);
        t_after = t(ch_t:end,:) - t(ch_t);
        y_after = lsim( sys_ch(ob_y_p,noise_point), u_after, t_after, x);
        y = [y_before;y_after(2:end,:)];
        if nargout > 1
            xhat_after = lsim( sys_ch(strcat('xhat_controlled',num2str(number)), noise_point), u_after, t_after, x);
            xhat = [xhat_before;xhat_after(2:end,:)];
        end
        if nargout > 2
            v_before = lsim( sys_ori(strcat('v_node',num2str(control_nodes1{number})),noise_point),u(1:ch_t,:),t(1:ch_t));
            v_after = lsim( sys_ch(strcat('v_node',num2str(control_nodes1{number})),noise_point), u_after, t_after, x);
            v = [v_before;v_after(2:end,:);];
        end
        if nargout > 3
            vhat_con_before = lsim( sys_ori(strcat('vhat_controlled',num2str(number)),noise_point),u(1:ch_t,:),t(1:ch_t));
            vhat_con_after = lsim( sys_ch(strcat('vhat_controlled',num2str(number)),noise_point), u_after, t_after, x);
            vhat_con = [vhat_con_before;vhat_con_after(2:end,:);];
        end
        if nargout > 4
            w_before = lsim( sys_ori(strcat('w_node',num2str(control_nodes1{number})),noise_point),u(1:ch_t,:),t(1:ch_t));
            w_after = lsim( sys_ch(strcat('w_node',num2str(control_nodes1{number})),noise_point), u_after, t_after, x);
            w = [w_before;w_after(2:end,:);];
        end
    end
    % Remove Controller
    if (Node1 == Node2) && (controller1 > controller2)
        % controller number is not remove number
        if number ~= condition
            % Add Controller
            y_before = lsim(sys_ori(ob_y_p,noise_point),u(1:ch_t,:),t(1:ch_t));
            x = lsim(sys_ori({'x'},noise_point),u(1:ch_t,:),t(1:ch_t));
            % Remove Controller State
            key = 1:controller1;
            key(key == condition) = [];
            for itr = key
                % Extend and Simple
                x_o_p = lsim(sys_ori(strcat('xhat_controlled',num2str(itr)), noise_point), u(1:ch_t,:), t(1:ch_t));
                x = [x(end,:),x_o_p(end,:)];
                if itr == number && nargout > 1
                    xhat_before = x_o_p;
                end
                % only Extend
                if controller_pattern1(itr) == 'e'
                    x_c_s = lsim(sys_ori( strcat('controller_x',num2str(itr)), noise_point), u(1:ch_t,:), t(1:ch_t));
                    x = [x,x_c_s(end,:)];
                end
            end 
            u_after = u(ch_t:end,:);
            t_after = t(ch_t:end,:) - t(ch_t);
            y_after = lsim( sys_ch(ob_y_p,noise_point), u_after, t_after, x);
            y = [y_before;y_after(2:end,:)];
            if nargout > 1
                xhat_after = lsim( sys_ch(strcat('xhat_controlled',num2str(number)), noise_point), u_after, t_after, x);
                xhat = [xhat_before;xhat_after(2:end,:)];
            end
            if nargout > 2
                v_before = lsim( sys_ori(strcat('v_node',num2str(control_nodes1{number})),noise_point),u(1:ch_t,:),t(1:ch_t));
                v_after = lsim( sys_ch(strcat('v_node',num2str(control_nodes1{number})),noise_point), u_after, t_after, x);
                v = [v_before;v_after(2:end,:);];
            end
            if nargout > 3
                vhat_con_before = lsim ( sys_ori(strcat('vhat_controlled',num2str(number)),noise_point),u(1:ch_t,:),t(1:ch_t));
                vhat_con_after = lsim( sys_ch(strcat('vhat_controlled',num2str(number)),noise_point), u_after, t_after, x);
                vhat_con = [vhat_con_before;vhat_con_after(2:end,:);];
            end
            if nargout > 4
                w_before = lsim( sys_ori(strcat('w_node',num2str(control_nodes1{number})),noise_point),u(1:ch_t,:),t(1:ch_t));
                w_after = lsim( sys_ch(strcat('w_node',num2str(control_nodes1{number})),noise_point), u_after, t_after, x);
                w = [w_before;w_after(2:end,:);];
            end
        else
            disp('Number_th controller is removed')
        end
    end
    % Add edge & wermove edge
    if (Node1 == Node2) && (controller1 == controller2)
        [ y_before, ~, x] = lsim(sys_ori(ob_y_p,noise_point),u(1:ch_t,:),t(1:ch_t));
        x = x(end, :);
        u_after = u(ch_t:end,:);
        t_after = t(ch_t:end,:) - t(ch_t);
        y_after = lsim( sys_ch(ob_y_p,noise_point), u_after, t_after, x);
        y = [y_before;y_after(2:end,:)];
        if nargout > 1
            xhat_before = lsim( sys_ori(strcat('xhat_controlled',num2str(itr)), noise_point), u(1:ch_t,:), t(1:ch_t));
            xhat_after = lsim( sys_ch(strcat('xhat_controlled',num2str(number)), noise_point), u_after, t_after, x);
            xhat = [xhat_before;xhat_after(2:end,:)];
        end
        if nargout > 2
            v_before = lsim( sys_ori(strcat('v_node',num2str(control_nodes1{number})),noise_point),u(1:ch_t,:),t(1:ch_t));
            v_after = lsim( sys_ch(strcat('v_node',num2str(control_nodes1{number})),noise_point), u_after, t_after, x);
            v = [v_before;v_after(2:end,:);];
        end
        if nargout > 3
            vhat_con_before = lsim( sys_ori(strcat('vhat_controlled',num2str(number)),noise_point),u(1:ch_t,:),t(1:ch_t));
            vhat_con_after = lsim( sys_ch(strcat('vhat_controlled',num2str(number)),noise_point), u_after, t_after, x);
            vhat_con = [vhat_con_before;vhat_con_after(2:end,:);];
        end
        if nargout > 4
            w_before = lsim( sys_ori(strcat('w_node',num2str(control_nodes1{number})),noise_point),u(1:ch_t,:),t(1:ch_t));
            w_after = lsim( sys_ch(strcat('w_node',num2str(control_nodes1{number})),noise_point), u_after, t_after, x);
            w = [w_before;w_after(2:end,:);];
        end
    end
end