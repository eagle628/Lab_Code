classdef Retrofit
    %RETROFIT_CONTROLLER 
    
    properties(Constant)
        Q_other = 0
    end
    
    methods
        
    end
    
    methods(Static)
        
        % Generate preexisting system 
        function sys = generate_fb(sys_local, sys_env)
            % generate_fb_Q(sys_local, Q)
            
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
            Bnew = [sys_local.b; zeros(m, size(sys_local.b, 2))];
            Cnew = [[C; G; DE*G], [tools.zeros([C;G], CE); CE]];
            sys = ss(Anew, Bnew, Cnew, 0);
            
            sys.InputGroup = sys_local.InputGroup;
            sys.OutputGroup.y = 1:size(C, 1);
            sys.OutputGroup.w = size(C, 1)+(1:size(G, 1));
            sys.OutputGroup.v = size(C, 1)+size(G,1)+(1:size(CE, 1));
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
            Bnew = [B, eye(size(B, 1), n); zeros(m, size(B, 2)+n)];
            Cnew = [[C; G; DE*G], [tools.zeros([C;G], CE); CE]; eye(size(A, 1), size(Anew, 1))];
            sys = ss(Anew, Bnew, Cnew, 0);
            sys.InputGroup.u = 1:size(B, 2);
            sys.InputGroup.d = size(B, 2)+ (1:n);
            sys.OutputGroup.y = 1:size(C, 1);
            sys.OutputGroup.w = size(C, 1)+(1:size(G, 1));
            sys.OutputGroup.v = size(C, 1)+size(G,1)+(1:size(CE, 1));
            sys.OutputGroup.z = size(C, 1)+size(G,1)+size(CE, 1)+(1:n);
        end
        

        
        function [sys, sys_K, gam] = design_hinf(sys_local, sys_model)
            sys_design = Retrofit.generate_fb_P(sys_local, sys_model);
            nm = size(sys_design({'y', 'w', 'v'}, :).c, 1);
            nc = size(sys_design(:, 'u').b, 2);
            [sys_K, sys_fb, gam] = hinfsyn(sys_design({'z', 'y', 'w', 'v'}, {'d', 'u'}), nm, nc, 'display', 'off');
            gam = norm(sys_fb, inf);
            sys = Retrofit.generate_retrofit_controller_dynamic(sys_local, sys_model, sys_K);
        end
        
        function [sys, sys_K] = design_h2(sys_local, sys_model)
            sys_design = Retrofit.generate_fb_P(sys_local, sys_model);
            nm = size(sys_design({'y', 'w', 'v'}, :).c, 1);
            nc = size(sys_design(:, 'u').b, 2);
            sys_K = h2syn(sys_design({'z', 'y', 'w', 'v'}, {'d', 'u'}), nm, nc);
            sys = Retrofit.generate_retrofit_controller_dynamic(sys_local, sys_model, sys_K);
        end
        
        % Use on test code in add.controller
        function [sys, K, sys_design, sys_design_K] = design(sys_local, sys_model, Q, R)
            if nargin < 4
                % Case : Not assign sys_model
                R = Q;
                Q = sys_model;
                % Design LQR
                K = lqr(sys_local('y', 'u'), Q, R);
                sys = Retrofit.generate_retrofit_controller_org(sys_local, K);
                sys_design = sys_local;
            else
                sys_design = Retrofit.generate_fb(sys_local, sys_model);
                Q1 = eye(size(sys_design.a))*Retrofit.Q_other;
                Q1(1:size(Q), 1:size(Q)) = Q;
                K = lqr(sys_design('y','u'), Q1, R);
                %                 sys_K = ss(sys_design.a-sys_design('y', 'u').b*K, sys_design.b, sys_design.c, sys_design.d);
                sys = Retrofit.generate_retrofit_controller_extend(sys_local, sys_model, K);
            end
            sys_design_K = Retrofit.connect_K(sys_design, K);
        end
        
        
        function sys = generate_retrofit_controller_dynamic(sys_local, sys_model, sys_K)
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
            B3 = Bk(:, ny+nw+(1:nv));
            D1 = Dk(:, 1:ny);
            D2 = Dk(:, ny+(1:nw));
            D3 = Dk(:, ny+nw+(1:nv));
            
            AA = [A+L*DE*G, -L*CE, tools.zeros(A, Ak);
                -BE*G, AE, tools.zeros(AE, Ak);
                -(B1*C+B2*G+B3*DE*G), B3*CE, Ak];
            BB = [tools.zeros(L, B1), -L*DE, L; tools.zeros(BE, B1), BE, tools.zeros(BE, L);
                B1, B2, tools.zeros(B2, L)];
            CC = [-D1*C-D2*G-D3*DE*G, D3*CE, Ck];
            DD = [D1, D2+D3*DE, tools.zeros(D1, L)];
            sys = ss(AA, BB, CC, DD);
            sys.InputGroup.y = 1:size(C, 1);
            sys.InputGroup.w = size(C, 1)+(1:size(G, 1));
            sys.InputGroup.v = size(C, 1)+size(G,1)+(1:size(CE, 1));
            sys.OutputGroup.u = 1:size(DD, 1);
        end
        
        % Extend Retro fit
        function sys = generate_retrofit_controller_extend(sys_local, sys_env, K)
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
            Bnew = [zeros(n+m, ny), [-L*DE, L; BE, tools.zeros(BE, L)]];
            Cnew = [K1, -K2; DE*G, CE; eye(size(A)), tools.zeros(A, K2)];
%             Cnew = [K1, -K2; 0*G, CE; eye(size(A)), tools.zeros(A, K2)];
            Dnew = [-K1, zeros(size(K1,1), nw+nv); zeros(size(K1, 1), size(Bnew, 2)); tools.zeros(A, Bnew)];
            
            sys = ss(Anew, Bnew, Cnew, Dnew);
            sys.InputGroup.y = 1:size(C, 1);
            sys.InputGroup.w = size(C, 1)+(1:size(G, 1));
            sys.InputGroup.v = size(C, 1)+size(G,1)+(1:size(CE, 1));
            sys.OutputGroup.u = 1:size(K, 1);
            sys.OutputGroup.vhat = size(K, 1) + (1:size(CE, 1));
            sys.OutputGroup.xhat = size(K, 1) + size(CE, 1) + (1:size(A, 1));
        end
        % Retro fit
        function sys = generate_retrofit_controller_org(sys_local, K)
            A = sys_local.a;
            Bv = sys_local(:,'v').b;    % \Gamma
            Cw = sys_local('w', :).c;   % Don't use
            Cy = sys_local('y', :).c;
            
            ny = size(Cy, 1);
            n = size(A, 1);
            nv = size(Bv, 2);
            
            % Ishizaki's paper (21) (output(フィードバックかける前提だからマイナスしている？))
            % C matrix of original system = identity
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

