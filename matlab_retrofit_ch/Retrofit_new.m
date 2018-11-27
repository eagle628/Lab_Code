classdef Retrofit_new
    %RETROFIT_CONTROLLER こ�??��クラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Constant)
        Q_other = 0
    end
    
    methods
        
    end
    
    methods(Static)
        
        function [sys, sys_design_K] = connect_controller(sys_all, sys_controller, sys_model,  Q, R)
            if nargin==5
                sys_local = sys_controller;
                [sys_controller, ~, ~, ~, sys_design_K] = Retrofit_new.design(sys_local, sys_model, Q, R);
                
            elseif nargin==4
                R = Q;
                Q = sys_model;
                sys_local = sys_controller;
                [sys_controller, ~, ~, ~, sys_design_K] = Retrofit_new.design(sys_local, Q, R);
            else
                sys_design_K = [];
            end
            
            variables = fieldnames(sys_controller.InputGroup);
            sys_controller = sys_controller(:, variables);
            A1 = sys_all.a;
            B1 = sys_all(:, 'u').b;
            C1 = sys_all(variables, :).c;
            
            A2 = sys_controller.a;
            B2 = sys_controller.b;
            C2 = sys_controller('u', :).c;
            D2 = sys_controller('u', :).d;
            
            Anew = [A1+B1*D2*C1, B1*C2; B2*C1, A2];
            Bnew = [sys_all.b; zeros(size(A2, 1), size(sys_all.b, 2))];
            Cnew = [sys_all.c, zeros(size(sys_all.c, 1), size(A2, 1));
                sys_controller(:,variables).d* C1, sys_controller(:,variables).c];
            sys = ss(Anew, Bnew, Cnew, 0);
            sys.InputGroup = sys_all.InputGroup;
            sys.OutputGroup = sys_all.OutputGroup;
            idx = size(sys_all.c, 1);
            keys = fieldnames(sys_controller.OutputGroup);
            for itr = 1:numel(keys)
                key = keys{itr};
                sys.OutputGroup.(key) = sys_controller.OutputGroup.(key) + idx;
                %                 idx = idx + numel(sys_controller.OutputGroup.(key));
            end
        end
        
        function sys = generate_fb(sys_local, sys_env)
            sys = Retrofit_new.generate_fb_P(sys_local, sys_env);
        end
        
        function sys = generate_fb_P(sys_local, sys_env)
            A = sys_local.a;
            L = sys_local(:,'v').b;
            G = sys_local('w', :).c;
            B = sys_local(:,'u').b;
            C = sys_local('y', :).c;
            AE= sys_env.A;
            BE = sys_env.B;
            CE = sys_env.C;
            DE = sys_env.D;
            n = size(A,1);
            m = size(AE,1);
            Anew = [A+L*DE*G, L*CE; BE*G, AE];
            
            if isfield(sys_local.InputGroup, 'd')
                Bd = sys_local(:, 'd').b;
            else
                Bd = eye(size(B, 1), n);
            end
            
            Bnew = [B, Bd; zeros(m, size(B, 2)+size(Bd, 2))];
            if isfield(sys_local.OutputGroup, 'z')
                Cz = sys_local('z', :).c;
                Cz = [Cz, tools.zeros(Cz, CE)];
            else
                Cz = eye(size(A, 1), size(Anew, 1));
            end
            
            if isfield(sys_local.OutputGroup, 'z') && isfield(sys_local.InputGroup, 'd')
                Dd = sys_local('z', 'd').d;
                Dnew = blkdiag(tools.zeros([C; G; DE*G], B), Dd);
            else
                Dnew = 0;
            end
            Cnew = [[C; G; DE*G], [myzeros([C;G], CE); CE]; Cz];
            sys = ss(Anew, Bnew, Cnew, Dnew);
            sys.InputGroup.u = 1:size(B, 2);
            sys.InputGroup.d = size(B, 2)+ (1:size(Bd, 2));
            sys.OutputGroup.y = 1:size(C, 1);
            sys.OutputGroup.w = size(C, 1)+(1:size(G, 1));
            sys.OutputGroup.v = size(C, 1)+size(G,1)+(1:size(CE, 1));
            sys.OutputGroup.z = size(C, 1)+size(G,1)+size(CE, 1)+(1:size(Cz, 1));
        end
        
        
        
        function [sys, sys_K, sys_design, sys_fb, gam] = design_dynamic(sys_local, sys_model, ratio, type)
            if nargin < 2 || isempty(sys_model)
                sys_model = ss(zeros(fliplr(size(sys_local('w', 'v')))));
            end
            if nargin < 3
                ratio = 0;
            end
            if nargin < 4
                type = 'H2';
            end
            sys_design = Retrofit_new.generate_fb_P(sys_local, sys_model);
            
            nc = size(sys_design(:, 'u').b, 2);
            direct = [zeros(nc, size(sys_design(:, 'd'), 2)), ratio*eye(nc)];
            
            C = sys_design({'y', 'w'}, :).c;
            D = sys_design({'y', 'w'}, {'d', 'u'}).d;
            CD = [C, D];
            
            [~, ia] = unique(CD, 'rows');
            idx_y = [sys_design.OutputGroup.y, sys_design.OutputGroup.w];
            idx_y = idx_y(ia);
            sys_design.OutputGroup.y_design = idx_y;
            nm = size(sys_design({'y_design'}, :).c, 1);
            if ~isinf(ratio)
                if ratio ~= 0
                    gplant = [direct; sys_design({'z', 'y_design'}, {'d', 'u'})];
                    gplant.OutputGroup.zu = 1;
                else
                    gplant = sys_design({'z', 'y_design'}, {'d', 'u'});
                end
                if strcmp(type, 'H2')
                    [sys_K, sys_fb, gam, info] = h2syn(gplant, nm, nc);
                else
                    %             [sys_K, sys_fb, gam, info] = hinfsyn(gplant, nm, nc, 'display', 'on', 'method', 'ric');
                    [sys_K, sys_fb, gam, info] = hinfsyn(gplant, nm, nc, 'display', 'on', 'method', 'lmi');
                end
                %             [sys_K, sys_fb, gam, info] = hinfsyn(gplant, nm, nc, 'display', 'on', 'method', 'maxe');
                
                %             gam = norm(sys_fb, inf);
                Bnew = zeros(order(sys_K), size(C, 1));
                Bnew(:, idx_y) = sys_K.b;
                Dnew = zeros(size(sys_K.d, 1), size(C, 1));
                Dnew(:, idx_y) = sys_K.d;
                sys_K = ss(sys_K.a, Bnew, sys_K.c, Dnew);
            else
                sys_K = ss(zeros(nc, size(C, 1)));
                gam = [];
                sys_fb = sys_design;
            end
            [sys, sys_K] = Retrofit_new.generate_retrofit_controller_dynamic(sys_local, sys_model, sys_K);
        end
        
        function [sys, sys_K, sys_design, sys_fb, gam] = design_dynamic_servo(sys_local, sys_model, ratio, type)
            if nargin < 2 || isempty(sys_model)
                sys_model = ss(zeros(fliplr(size(sys_local('w', 'v')))));
            end
            if nargin < 3
                ratio = 0;
            end
            if nargin < 4
                type = 'H2';
            end
            sys_design = Retrofit_new.generate_fb_P(sys_local, sys_model);
            
            nc = size(sys_design(:, 'u').b, 2);
            

            
            sys_design = [zeros(size(sys_design, 1), order(sys_local)), sys_design;
                eye(order(sys_local)), zeros(order(sys_local), size(sys_design, 2))];
            sys_design.InputGroup.dx = 1:order(sys_local);
            sys_design.OutputGroup.dx = size(sys_design, 1)-order(sys_local)+1:size(sys_design, 1);
            sys_design('z', 'dx') = sys_design('z', 'dx') + sys_local('z', :).c * 1/(tf('s')+0.01);
            C = sys_design({'y', 'w', 'dx'}, :).c;
            D = sys_design({'y', 'w', 'dx'}, {'dx', 'd', 'u'}).d;
            CD = [C, D];
            [~, ia] = unique(CD, 'stable', 'rows');
            idx_y = 1:size([sys_design.OutputGroup.y, sys_design.OutputGroup.w, sys_design.OutputGroup.dx], 2);
            idx_y = idx_y(ia);
            sys_design.OutputGroup.y_design = idx_y;
            nm = size(sys_design({'y_design'}, :).c, 1);
            if ~isinf(ratio)
                if ratio ~= 0
                    direct = [zeros(nc, size(sys_design(:, {'dx', 'd'}), 2)), ratio*eye(nc)];
                    gplant = [direct; sys_design({'z', 'y_design'}, {'dx', 'd', 'u'})];
                    gplant.OutputGroup.zu = 1;
                else
                    gplant = sys_design({'z', 'y_design'}, {'d', 'u'});
                end
                if strcmp(type, 'H2')
                    [sys_K, sys_fb, gam, info] = h2syn(gplant, nm, nc);
                else
                    %             [sys_K, sys_fb, gam, info] = hinfsyn(gplant, nm, nc, 'display', 'on', 'method', 'ric');
                    [sys_K, sys_fb, gam, info] = hinfsyn(gplant, nm, nc, 'display', 'on', 'method', 'lmi');
                end
                %             [sys_K, sys_fb, gam, info] = hinfsyn(gplant, nm, nc, 'display', 'on', 'method', 'maxe');
                
                %             gam = norm(sys_fb, inf);
                Bnew = zeros(order(sys_K), size(C, 1));
                Bnew(:, idx_y) = sys_K.b;
                Dnew = zeros(size(sys_K.d, 1), size(C, 1));
                Dnew(:, idx_y) = sys_K.d;
                sys_K = ss(sys_K.a, Bnew, sys_K.c, Dnew);
            else
                sys_K = ss(zeros(nc, size(C, 1)));
                gam = [];
                sys_fb = sys_design;
            end
            [sys, sys_K] = Retrofit_new.generate_retrofit_controller_dynamic(sys_local, sys_model, sys_K);
        end
        
        
        function [sys, sys_K, sys_design, sys_fb, gam] = design_hinf2(sys_local, sys_model, ratio, ratio2)
            if nargin < 3
                [sys, sys_K, sys_design, sys_fb, gam] = design_hinf(sys_local, sys_model);
            else
                sys_design = Retrofit_new.generate_fb_P(sys_local, sys_model);
                Cnew = [sys_design('z', :).c; ratio*sys_design('v', :).c];
                Dnew = [sys_design('z', :).d; ratio*sys_design('v', :).d];
                if nargin == 4
                    Cnew = [Cnew; ratio2*sys_design('w', :).c];
                    Dnew = [Dnew; ratio2*sys_design('w', :).d];
                end
                sys_design2 = ss(sys_design.a, sys_design.b, [sys_design.c; Cnew], [sys_design.d; Dnew]);
                sys_design2.InputGroup = sys_design.InputGroup;
                sys_design2.OutputGroup = sys_design.OutputGroup;
                sys_design2.OutputGroup.z2 =size(sys_design.c, 1)+(1:size(Cnew, 1));
                nm = size(sys_design2({'y', 'w'}, :).c, 1);
                nc = size(sys_design2(:, 'u').b, 2);
                [sys_K, sys_fb, gam] = hinfsyn(sys_design2({'z2', 'y', 'w'}, {'d', 'u'}), nm, nc, 'display', 'off');
                %             [sys_K, sys_fb, gam] = h2syn(sys_design2({'z2', 'y', 'w'}, {'d', 'u'}), nm, nc);
                
                %             gam = norm(sys_fb, inf);
                [sys, sys_K] = Retrofit_new.generate_retrofit_controller_dynamic(sys_local, sys_model, sys_K);
            end
        end
        
        %         function [sys, sys_K] = design_h2(sys_local, sys_model)
        %             sys_design = Retrofit_new.generate_fb_P(sys_local, sys_model);
        %             nm = size(sys_design({'y', 'w', 'v'}, :).c, 1);
        %             nc = size(sys_design(:, 'u').b, 2);
        %             sys_K = h2syn(sys_design({'z', 'y', 'w', 'v'}, {'d', 'u'}), nm, nc);
        %             sys = Retrofit_new.generate_retrofit_controller_dynamic(sys_local, sys_model, sys_K);
        %         end
        
        function [sys, K, sys_design, sys_design_K] = design(sys_local, sys_model, Q, R, Qv, Qw, QxE)
            if nargin < 4
                R = Q;
                Q = sys_model;
                K = lqr(sys_local('y', 'u'), Q, R);
                sys = Retrofit_new.generate_retrofit_controller_org(sys_local, K);
                sys_design = sys_local;
            else
                if isempty(sys_model)
                    sys_model = ss(tf(zeros(fliplr(size(sys_local('w', 'v'))))));
                end
                sys_design = Retrofit_new.generate_fb(sys_local, sys_model);
                Q1 = eye(size(sys_design.a))*Retrofit_new.Q_other;
                Q1(1:size(Q), 1:size(Q)) = Q;
                if nargin > 4
                    Cv = sys_design('v',:).c;
                    Q1 = Q1 + Cv'*Qv*Cv;
                end
                if nargin > 5
                    Cw = sys_design('w', :).c;
                    Q1 = Q1 + Cw'*Qw*Cw;
                end
                if nargin > 6
                    Q2 = zeros(size(sys_design.a));
                    Q2(size(Q, 1)+1:end, size(Q,1)+1:end) = QxE;
                    Q1 = Q1 + Q2;
                end
                if isinf(R)
                    K = zeros(size(sys_design('y', 'u')))';
                else
                    K = lqr(sys_design(:,'u'), Q1, R);
                end
                %                 sys_K = ss(sys_design.a-sys_design('y', 'u').b*K, sys_design.b, sys_design.c, sys_design.d);
                sys = Retrofit_new.generate_retrofit_controller(sys_local, sys_model, K);
            end
            sys_design_K = Retrofit_new.connect_K(sys_design, K);
        end
        
        
        function [sys, sys_K] = generate_retrofit_controller_dynamic(sys_local, sys_model, sys_K)
            A = sys_local.a;
            L = sys_local(:,'v').b;
            G = sys_local('w', :).c;
            %             B = sys_local(:,'u').b;
            C = sys_local('y', :).c;
            
            ny = size(C, 1);
            nw = size(G, 1);
            nv = size(L, 2);
            
            AE= sys_model.A;
            BE = sys_model.B;
            CE = sys_model.C;
            DE = sys_model.D;
            
            [Ak ,Bk, Ck, Dk] = ssdata(sys_K);
            B1 = Bk(:, 1:ny);
            B2 = Bk(:, ny+(1:nw));
            % %             B3 = Bk(:, ny+nw+(1:nv));
            %             D1 = Dk(:, 1:ny);
            %             D2 = Dk(:, ny+(1:nw));
            %             D3 = Dk(:, ny+nw+(1:nv));
            
            AA = [A+L*DE*G, -L*CE, myzeros(A, Ak);
                -BE*G, AE, myzeros(AE, Ak);
                -Bk*[C;G], myzeros(Bk, AE), Ak];
            BB = [myzeros(L, B1), -L*DE, L;
                myzeros(BE, B1), BE, myzeros(BE, L);
                Bk, myzeros(B2, L)];
            CC = [-Dk*[C;G], myzeros(Dk, CE), Ck;
                DE*G, CE, myzeros(DE, Ck);
                eye(size(A, 1)), myzeros(A, CE), myzeros(A, Ck)];
            DD = [Dk, myzeros(Dk, L)];
            DD = [DD; myzeros(DE, DD); myzeros(A, DD)];
            sys = ss(AA, BB, CC, DD);
            sys.InputGroup.y = 1:size(C, 1);
            sys.InputGroup.w = size(C, 1)+(1:size(G, 1));
            sys.InputGroup.v = size(C, 1)+size(G,1)+(1:size(CE, 1));
            sys.OutputGroup.u = 1:size(Ck, 1);
            sys.OutputGroup.vhat = size(Ck, 1) + (1:size(CE, 1));
            sys.OutputGroup.xhat = size(Ck, 1) + size(CE, 1) + (1:size(A, 1));
            
            sys_K.OutputGroup.u = 1:size(Ck, 1);
            sys_K.InputGroup.y = 1:size(C, 1);
            sys_K.InputGroup.w = size(C, 1)+(1:size(G, 1));
        end
        
        function sys = generate_retrofit_controller(sys_local, sys_env, K)
            A = sys_local.a;
            L = sys_local(:,'v').b;
            G = sys_local('w', :).c;
            B = sys_local(:,'u').b;
            C = sys_local('y', :).c;
            
            AE= sys_env.A;
            BE = sys_env.B;
            CE = sys_env.C;
            DE = sys_env.D;
            
            ny = size(C, 1);
            n = size(A, 1);
            nw = size(G, 1);
            nv = size(L, 2);
            
            m = size(AE, 1);
            
            K1 = K(:, 1:ny);
            K2 = K(:, ny+1:end);
            %
            Anew = [A+L*DE*G, -L*CE; -BE*G, AE];
            Bnew = [zeros(n+m, ny), [-L*DE, L; BE, myzeros(BE, L)]];
            Cnew = [K1, -K2; DE*G, CE; eye(size(A)), myzeros(A, K2)];
            Dnew = [-K1, zeros(size(K1,1), nw+nv); myzeros(DE, Bnew); myzeros(A, Bnew)];
            
            sys = ss(Anew, Bnew, Cnew, Dnew);
            sys.InputGroup.y = 1:size(C, 1);
            sys.InputGroup.w = size(C, 1)+(1:size(G, 1));
            sys.InputGroup.v = size(C, 1)+size(G,1)+(1:size(CE, 1));
            sys.OutputGroup.u = 1:size(K, 1);
            sys.OutputGroup.vhat = size(K, 1) + (1:size(CE, 1));
            sys.OutputGroup.xhat = size(K, 1) + size(CE, 1) + (1:size(A, 1));
        end
        
        function sys = generate_retrofit_controller_org(sys_local, K)
            A = sys_local.a;
            Bv = sys_local(:,'v').b;
            Cw = sys_local('w', :).c;
            Cy = sys_local('y', :).c;
            
            ny = size(Cy, 1);
            n = size(A, 1);
            nv = size(Bv, 2);
            
            Anew = A;
            Bnew = [zeros(n, ny), Bv];
            Cnew = K;
            Dnew = [-K, zeros(size(K,1), nv)];
            
            sys = ss(Anew, Bnew, Cnew, Dnew);
            sys.InputGroup.y = 1:size(Cy, 1);
            sys.InputGroup.v = size(Cy, 1)+(1:nv);
            sys.OutputGroup.u = 1:size(Dnew, 1);
            
        end
        
        function sys_K = connect_K(sys_design, K)
            sys_K = ss(sys_design.a-sys_design('y', 'u').b*K, sys_design.b, [sys_design.c; -K], [sys_design.d; zeros(size(K, 1), size(sys_design.d, 2))]);
            sys_K.InputGroup = sys_design.InputGroup;
            sys_K.OutputGroup = sys_design.OutputGroup;
            sys_K.OutputGroup.u = size(sys_design.c, 1) + (1:size(K, 1));
        end
        %
        
    end
    
end
%
function z = myzeros(A, B)
z = zeros(size(A, 1), size(B, 2));
end
