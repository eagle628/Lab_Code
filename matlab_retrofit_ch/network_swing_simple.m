classdef network_swing_simple < handle
    %{
    Adj : Adjacency matrix (Finite Graph representation)
    Adj_ref : Wall Connect
    N : number of  Node
    seed : Randam number Seed
    nodes : All node Objects
    controllers :
    %}
    
    properties
        Adj
        Adj_ref
        N
        seed
        nodes
        controllers
    end
    
    methods
        % Generate Network
        % Y Garph Weight Control
        % The bigger "rc" becomes, the more graph is density.
        function obj = network_swing_simple(N, m, d, b, Y, rc, seed)
            % Nargin Missing
            if nargin < 2
                m = [];
            end
            if nargin < 3
                d = [];
            end
            if nargin < 4
                b = [];
            end
            if nargin < 5
                Y = [];
            end
            if nargin < 6
                rc = [];
            end
            r = rng();
            if nargin==0 || isempty(N)
                N = 30;
            end
            % Set Randam number
            if nargin == 7
                rng(seed);
                obj.seed = seed;
            else
                obj.seed = r;
            end
            % Set Node number
            obj.N = N;
            % Get Nodes Memory for object
            obj.nodes = cell(N, 1);
            % Get Nodes object
            for itr = 1:N
                obj.nodes{itr} = power_swing_2nd(m, d, b, strcat('g', num2str(itr)));
                for jj = 1:5, rand(); end
            end
            % Generate Network Structure
            obj.generate_graph(Y, rc, obj.seed);
        end
        
        % plot Network Structure
        function plot(obj)
            G = graph(obj.Adj);
            fnet = tools.myplot([], G, 'EdgeLabel',G.Edges.Weight)
            fnet.Name = 'Network';
        end
        
        % Add Input u&d, Output y,v&w to the idx_th subsystem 
        % idx が 1つのみしか受け入れないとする
        function sys_out = add_io(obj, sys, idx, name)
            if iscell(idx)
                % idx がセルならば
                sys_out = sys;
                for itr = 1:numel(idx)
                    idx_i = idx{itr};
                    if nargin == 4
                        name_i = name{itr};
                    else
                        name_i = num2str(itr);
                    end
                    sys_out = obj.add_io(sys_out, idx_i, name_i);
                end                
            else
                if nargin < 4
                    name = '';
                end
                if ~isempty(name)
                    name = strcat('_', name);
                end
                G_Lap = obj.get_Laplacian_out(idx);
                B1 = sys(:, 'u').b(:, idx);
                Cv = -G_Lap(idx, :)*sys('theta',:).c;
                Cw = sys('theta',:).c(idx, :);
                IDX = [idx(:)*2-1, idx(:)*2]';
                Cy = sys('x',:).c(IDX(:),:);
                Bd = sys(:, 'x').b(:, IDX(:));
                sys_out = ss(sys.a, [sys.b, B1, Bd], [sys.c; Cy; Cv; Cw], blkdiag(sys.d, zeros(size([Cy; Cv; Cw], 1), size([B1, Bd], 2))));
                sys_out.InputGroup = sys.InputGroup;
                sys_out.OutputGroup = sys.OutputGroup;
                % add name input
                if isstring(name)
                    name = char(name);
                end
                % Control Input
                base = size(sys.b, 2); 
                sys_out.InputGroup.(strcat('u', name)) = base+(1:size(B1, 2));
                % Disturbance
                base = base + size(B1, 2);
                sys_out.InputGroup.(strcat('d', name)) = base+(1:size(Bd, 2));
                % add name output
                % Mesurement
                base = size(sys.c, 1);
                sys_out.OutputGroup.(strcat('y', name)) = base + (1:size(Cy, 1));
                % Interconnection Feed Back 
                base = base + size(Cy, 1);
                sys_out.OutputGroup.(strcat('v', name)) = base + (1:size(Cv, 1));
                % Interconnection Output
                base = base + size(Cv, 1);
                sys_out.OutputGroup.(strcat('w', name)) = base + (1:size(Cw, 1));
            end
            
        end
        
        % Remove_node
        % Keep Node Number
        % Idx is array.
        function remove_node(obj, idx)
%             nidx = obj.get_not_idx(idx);
%             obj.Adj = obj.Adj(nidx, nidx);
%             obj.nodes(idx) = [];
            Mat = obj.Adj;
            Mat( idx ,:) = 0;
            Mat( :, idx) = 0;
            obj.Adj = Mat;
%             obj.N = obj.N-numel(idx);
        end
        
        % Add_node
        function add_node(obj,m,d,b,Y,ratio_connect,seed,ref)
            if nargin < 8
                ref = 0;
            end
            if nargin < 7
                seed = rng();
                seed = seed.Seed;
            end
            % New Node NUmber is obj.N+1.
            % Make New Node prameter
            obj.nodes{obj.N+1} = power_swing_2nd(m,d,b, strcat('add', num2str(obj.N+1)),seed);
            
            % Add Adj matrix
            thre_edge = 1-ratio_connect;
            connect = 0;
            count = 0;
            % Except Isolation Node
            D = sum(obj.Adj);
            idx = find(D~=0);
            Adj_mat = obj.Adj(idx,idx);
            while ~connect
                count=count+1;
                add_adj = rand(length(Adj_mat),1);
                edge_weight = rand(length(Adj_mat),1);
                add_adj = double(add_adj>thre_edge);
                add_adj = (Y(1) + edge_weight*(Y(2)-Y(1))).*add_adj;
                new_Adj = [Adj_mat,add_adj;add_adj',0;];
                D_mat = diag(sum(new_Adj,2));
                Lap_G = D_mat-new_Adj;
                connect = sum(eig(Lap_G)<1e-5)==1;
            end
            % set new Network properties
            obj.Adj_ref = [obj.Adj_ref; ref];
            obj.N = obj.N + 1; 
            Adj_mat = zeros(obj.N);
            Adj_mat([idx,obj.N],[idx,obj.N]) = new_Adj;
            obj.Adj = Adj_mat;
        end
        
        % Add Edge Or Edge Weght OverRide
        function add_edge(obj,idx2,Y,seed)
            if nargin < 4
                rng(6);
            else
                rng(seed);
            end
            if nargin < 3
                Y = [1,5];
            end
            if numel(Y) == 1
                Y = [Y, Y];
            end
            % Except Isolation Node
            D = sum(obj.Adj);
            idx = find(D~=0);
            if isempty(idx2)
                idx2 = [idx(1),idx(2)];
            end
            Adj_mat = obj.Adj(idx,idx);
            connect = 0;
            count = 0;
            while ~connect
                count=count+1;
                edge_weight = rand(1);
                add_edge = Y(1) + edge_weight*(Y(2)-Y(1));
                Adj_mat(idx==idx2(1), idx==idx2(2)) = add_edge;
                Adj_mat(idx==idx2(2), idx==idx2(1)) = add_edge;
                D_mat = diag(sum(Adj_mat,2));
                Lap_G = D_mat-Adj_mat;
                connect = sum(eig(Lap_G)<1e-5)==1;
            end
            new_Adj = zeros(obj.N);
            new_Adj(idx,idx) = Adj_mat;
            obj.Adj = new_Adj;
        end
        
        % Remove Edge
        function remove_edge(obj,idx2,seed)
            if nargin < 3
                rng(6);
            else
                rng(seed);
            end
            obj.Adj(idx2(1),idx2(2)) = 0;
            obj.Adj(idx2(2),idx2(1)) = 0;
        end     
            
        % Set Grpaph Structre (Don't Use)
        function set_Adj(obj, i, j, Y)
            obj.Adj(i,j) = Y;
            obj.Adj(j, i) = Y;
        end
        
        % changem partial parameter for m
        function set_M(obj, i, Ms)
            if numel(i) ~= numel(Ms)
                warning('Invalid.\n');
            else
                for itr = 1:numel(i)
                    obj.nodes{itr}.m = Ms(itr);
                end
            end
        end
        
        % change partial parameter for d
        function set_D(obj, i, Ds)
            if numel(i) ~= numel(Ds)
                warning('Invalid.\n');
            else
                for itr = 1:numel(i)
                    obj.nodes{itr}.d = Ds(itr);
                end
            end
        end
        
        % Get partial parameter for m
        function Ms = get_M(obj, idx)
            if nargin < 2
                idx = 1:obj.N;
            end
            Ms = zeros(numel(idx), 1);
            for itr = 1:numel(idx)
                Ms(itr) = obj.nodes{idx(itr)}.m;
            end
            
        end
        
        % Get partial parameter for d
        function Ds = get_D(obj, idx)
            if nargin < 2
                idx = 1:obj.N;
            end
            Ds = zeros(numel(idx), 1);
            for itr = 1:numel(idx)
                Ds(itr) = obj.nodes{idx(itr)}.d;
            end
            
        end
        
        % Generate Laplacian Matrix
        function Lap_G = get_Laplacian(obj, idx, ref)
            % Nargine Missing
            if nargin < 2
                Adj_Mat = obj.Adj;
                Adj_r = obj.Adj_ref;
            else
                % Extract local system relation
                Adj_Mat = obj.Adj(idx, idx);
                if ref
                    Adj_r = obj.Adj_ref(idx);
                else
                    Adj_r = 0;
                end
            end
            D_mat = diag(sum(Adj_Mat,2)+Adj_r);
            Lap_G = D_mat-Adj_Mat;
        end
        
        % Remove Partition Laplacian Matrix
        function Lap_G = get_Laplacian_out(obj, idx)
            Adj_Mat = obj.Adj;
            for itr = 1:numel(idx)
                for jtr = itr:numel(idx)
                    Adj_Mat(idx(itr), idx(jtr)) = 0;
                    Adj_Mat(idx(jtr), idx(itr)) = 0;
                end
            end
            Adj_r = obj.Adj_ref;
            D_mat = diag(sum(Adj_Mat,2)+Adj_r);
            Lap_G = D_mat-Adj_Mat;
        end
        
        % Get controlled system
        % add_controller を行った後にしか稼働しない．
        % sys_all は add_ioなどで，出力ポートの付加を行った後のコントローラをつける段階のもの
        function sys_controlled = get_sys_controlled(obj, sys_all)
            for itr = 1:numel(obj.controllers)
                %idx = obj.names2idx(obj.controllers{itr}.names);
                idx = obj.controllers{itr}.nodes;
                G_Lap = obj.get_Laplacian_out(idx);
                A1 = sys_all.a;
                B1 = sys_all(:, 'u').b(:, idx);
                Cv = -G_Lap(idx, :)*sys_all('theta',:).c;   % Flip interconnection affect 
                Cw = sys_all('theta',:).c(idx, :);
                IDX = [idx(:)*2-1, idx(:)*2]';
                Cy = sys_all('x',:).c(IDX(:),:);
                C1 = [Cy; Cw; Cv];
                sys_controller = obj.controllers{itr}.sys;
                [A2, B2, C2, D2] = ssdata(sys_controller('u', {'y', 'w', 'v'}));
                Anew = [A1+B1*D2*C1, B1*C2; B2*C1, A2];
                Bnew = [sys_all.b; zeros(size(A2, 1), size(sys_all.b, 2))];
                Cnew = [sys_all.c, zeros(size(sys_all.c, 1), size(A2, 1));
                    sys_controller({'u','vhat'},{'y', 'w', 'v'}).d* C1, sys_controller({'u','vhat'},{'y', 'w', 'v'}).c;
                    zeros(length(A2),length(A1)),eye(length(A2))];
                sys = ss(Anew, Bnew, Cnew, 0);
                sys.InputGroup = sys_all.InputGroup;
                sys.OutputGroup = sys_all.OutputGroup;
                idx = size(sys_all.c, 1);
                keys = fieldnames(sys_controller.OutputGroup);
                % C_newにおいては，コントローラーの状態について一気につけているが，
                % output rectiferとextendによるコントローラの状態について名づけを別にしている．
                for itr2 = 1:numel(keys)
                    key = keys{itr2};
                    sys.OutputGroup.(strcat(key, '_controlled', num2str(itr))) = sys_controller.OutputGroup.(key) + idx;
                    %                 idx = idx + numel(sys_controller.OutputGroup.(key));
                end
                idx = idx + sys_controller.OutputGroup.(key)(end);
                % このlineがretro(Not Extend)で動くことは確認済み
                sys.OutputGroup.(strcat('controller_x',num2str(itr))) = idx+1:idx+length(A2)-numel(sys_controller.OutputGroup.(key)); 
                sys_all = sys;   
            end
            sys_controlled = sys_all;
            % Controller Infomation
            sys_controlled.Userdata.controllers = numel(obj.controllers);
            sys_controlled.Userdata.Nodes = obj.N;
            sys_controlled.Userdata.controllertype = '';
            sys_controlled.Userdata.control_nodes = {};
            for itr = 1:numel(obj.controllers)
                sys_controlled.Userdata.controllertype = ...
                            strcat(sys_controlled.Userdata.controllertype,obj.controllers{itr}.type);
                sys_controlled.Userdata.control_nodes{itr} = obj.controllers{itr}.nodes;
            end
        end
               
        % Add controller for idx_th subsystem
        % Extend Retro fit OR Retro fit is(model)
        function add_controller(obj, idx, model, Q, R)
            sys_local = obj.get_sys_local(idx);
            nu = size(sys_local(:, 'u').b, 2);
            flag = false;
            if nargin == 2
                model = ss([], [], [], sys_local('w', 'v').d*0);
                flag = true;
            end
            if nargin == 2 || nargin==3
                Q = kron(eye(numel(idx)), diag([1, 100]));
                R = eye(nu);
            elseif nargin == 4
                R = Q;
                Q = model;
                model = ss([], [], [], sys_local('w', 'v').d*0);
                flag = true;
            end
            if isempty(Q)
                Q = kron(eye(numel(idx)), diag([1, 100]));
            end
            if isempty(R)
                R = eye(nu);
            end
            if isempty(model)
                model = ss([], [], [], sys_local('w', 'v').d*0);
                flag = true;
            end
            
            [sys_control, K, sys_design, sys_K] = Retrofit.design(sys_local, model, Q, R);
            controller = struct();
            controller.sys = sys_control;
            controller.K = K;
            controller.sys_design = sys_design;
            controller.sys_K = sys_K;
            controller.Q = Q;
            controller.R = R;
            controller.model = model;
            controller.names = {};
            % controller type
            if flag
                controller.type = 's';
            else
                controller.type = 'e';
            end
            % controller name
            for itr = 1:numel(idx)
                controller.names = [controller.names; {obj.nodes{idx(itr)}.name}];
            end
            % control node
            controller.nodes = idx;
            obj.controllers = [obj.controllers; {controller}];
        end
        
        % ゲインを変えずにモデルだけ変更(Extend )
        function add_controller2(obj, idx, model, K)
            sys_local = obj.get_sys_local(idx);
            nu = size(sys_local(:, 'u').b, 2);
            
            [sys_control, K, sys_design, sys_K] = Retrofit.design2(sys_local, model, K);
            controller = struct();
            controller.sys = sys_control;
            controller.K = K;
            controller.sys_design = sys_design;
            controller.sys_K = sys_K;
            controller.model = model;
            controller.names = {};
            % controller type
            controller.type = 'e';
            % controller name
            for itr = 1:numel(idx)
                controller.names = [controller.names; {obj.nodes{idx(itr)}.name}];
            end
            % control node
            controller.nodes = idx;
            obj.controllers = [obj.controllers; {controller}];
        end
        
        % get_sys_controlled 内で，controller node名をとるためかな
        function idx = names2idx(obj, names)
            [c, ia, ic] = unique(names);
            idx = nan(numel(c), 1);
            for itr = 1:numel(obj.nodes)
                [f, i] = ismember(obj.nodes{itr}.name, c);
                if f
                    idx(i) = itr;
                end
            end
            idx = idx(ic);
        end
        
        % remove unnecessary node idx 
        function idx_out = get_not_idx(obj, idx)
            idx_out = (1:obj.N)';
            idx_out(idx) = [];
        end
        
        % Divide the system with idx to local system and environmental system
        function [sys_local, sys_env] = get_sys_local(obj, idx)
            n_local = numel(idx);
            % idxにおけるローカルシステムの取り出し
            sys = obj.get_sys(idx);
            % ローカルシステムにw,yを付け足す．
            A = sys.a;
            Cw = kron(eye(n_local), [1, 0]);
            Cy = eye(size(A));
            sys_local = ss(A, sys.b, [sys.c; Cy; Cw], 0);
            sys_local.InputGroup = sys.InputGroup;
            sys_local.OutputGroup = sys.OutputGroup;
            sys_local.OutputGroup.y = size(sys.c, 1) + (1:size(Cy, 1));
            sys_local.OutputGroup.w = size(sys.c, 1) + size(Cy,1) + (1:size(Cw, 1));

            if nargout == 2
                nidx = obj.get_not_idx(idx);
                con = obj.Adj(nidx, idx);
                J = -diag(sum(con)'+obj.Adj_ref(idx)); % interconnection signal -Y_{ij}*\theta_{i} in v 
                sys_e = obj.get_sys(obj.get_not_idx(idx), true);
                F = sys_e.a;
                if isempty(F)
                    G = [];
                    H = [];
                else
%                     Y = obj.Adj(obj.get_not_idx(idx), idx);

                    Y = obj.Adj(nidx, idx);
                    G = sys_e(:, 'v').b*Y; % extract input v term (summarize by n_local) (otherwise "w")
                    H = Y'*sys_e('theta', :).c; % interconnection signal Y_{ij}*\theta_{j} in v
                    D = -diag(sum(con, 2)); % J?
                    F = F + sys_e(:, 'v').b*D*sys_e('theta', :).c;
                end
                sys_env = ss(F, G, H, J);
            end
            
        end
        
        % Propertoes (Adjacency matrix)
        function generate_graph(obj, Y, ratio_connect, seed)
            % Generate Network Structure
            % Set Randam Numer
            s = rng();
            if nargin == 4
                rng(seed);
            end
            % Nargin Missing
            if nargin < 3 || isempty(ratio_connect)
                ratio_connect = 0.1;
            end
            if nargin < 2 || isempty(Y)
                Y = [1, 10];
            end
            
            if numel(Y) == 1
                Y = [Y, Y];
            end
            count=0;
            connect = false;
            N = obj.N; %#ok
            thre_edge = 1-ratio_connect;
            while ~connect
                count=count+1;
                Adj_seed = rand(N); %#ok
                Adj_seed2 = rand(N); %#ok
                Adj_ori = zeros(N); %#ok
                Adj_ori(Adj_seed>thre_edge) = 1;
                Adj_ori(Adj_seed<=thre_edge) = 0;
                Adj_ori = (Y(1) + Adj_seed2*(Y(2)-Y(1))).*Adj_ori;
                Adj_Mat = triu(Adj_ori,1);
                Adj_Mat = Adj_Mat+Adj_Mat';
                D_mat = diag(sum(Adj_Mat,2));
                Lap_G = D_mat-Adj_Mat;
                connect = sum(eig(Lap_G)<1e-5)==1;
%                 connect = sum(eig(Lap_G)<1e-10)==1;
%                 connect = (sum(eig(Lap_G)<1e-5) + sum(eig(Lap_G)<0))==1;
                    % Lap_G is semi-positive definit.
            end
            % Adjacency matrix
            obj.Adj = Adj_Mat;
            obj.Adj_ref = zeros(N, 1); %#ok
            obj.Adj_ref(end) = Y(1) + rand*(Y(2)-Y(1));
            rng(s);
        end
        
        % System in the area specified by "idx"
        function sys = get_sys(obj, idx, ref)
            % Nargin Missing
            if nargin < 2
                Lap_G = obj.get_Laplacian();
                idx = 1:obj.N;
            else
                if nargin < 3
                    ref = false;
                end
                Lap_G = obj.get_Laplacian(idx, ref);
            end
            % Get Memory
            Ms = zeros(numel(idx), 1);
            Ds = zeros(numel(idx), 1);
            for itr = 1:numel(idx)
                Ms(itr) = obj.nodes{idx(itr)}.m;
                Ds(itr) = obj.nodes{idx(itr)}.d;
            end
            % Node Parameter 
            Ms = diag(Ms);
            Ds = diag(Ds);
            
            % State Space A,B,C,D
            A = kron(eye(numel(idx)),[0 1;0 0])-kron(Ms\Ds,[0 0; 0 1])-kron(Ms\Lap_G,[0 0;1 0]);
            B = [];
            for itr = 1:numel(idx)
                B = blkdiag(B, [0; obj.nodes{idx(itr)}.b/obj.nodes{idx(itr)}.m]);
            end
            Bv = [];
            for itr = 1:numel(idx)
                Bv = blkdiag(Bv, [0; 1/obj.nodes{idx(itr)}.m]);
            end
            Bd = eye(size(A));
            %             if nargin < 2
            %                 Anew = A;
            %                 Anew(end, end-1) = Anew(end, end-1)-1;
            %                 A = Anew;
            %             end
            sys = ss(A, [B, Bv, Bd], eye(size(A)), 0);
            sys.OutputGroup.x = 1:2*numel(idx);
            sys.OutputGroup.omega = 2:2:2*numel(idx);
            sys.OutputGroup.theta = 1:2:2*numel(idx);
            sys.InputGroup.u = 1:numel(idx);
            sys.InputGroup.v = numel(idx)+(1:numel(idx));
            sys.InputGroup.x = numel(idx)*2+(1:numel(idx)*2);
            %             for itr = 1:numel(idx)
            %                 sys.InputGroup.(strcat('u_', obj.nodes{idx(itr)}.name)) = itr;
            %             end
            %             for itr = 1:numel(idx)
            %                 sys.OutputGroup.(strcat('y_', obj.nodes{idx(itr)}.name)) = itr*2-[1, 0];
            %             end
        end
        
        % ?
        function varargout = initial(obj, x0in, varargin)
            if nargin == 1 || isempty(x0in)
                x0in = [1; 1];
            end
            sys = obj.get_sys_controlled();
            varargout = cell(nargout, 1);
            x0 = zeros(size(sys.a, 1), 1);
            for itr = 1:numel(obj.nodes)
                x0(2*itr-[1, 0]) = x0in;
            end
            [varargout{:}] = initial(sys('omega', :), x0, varargin{:});
        end
        
        % model error
        % not complete
        function sys_frd = model_error(obj,sys,idx,e_r)
            sys_frd = sys;
            if iscell(idx)
                idx = cell2mat(idx);
            end
            A = sys_frd.A;
            k1 = ureal('k1',1,'Per',[-20 20]);
            trans_A = eye(size(A,1));
            for i = 1 : numel(idx)
                trans_A(idx(i),idx(i)) = usample(k1,1);
            end
            A=trans_A*A;
            sys_frd.A = A;
            
        end
    end
end



