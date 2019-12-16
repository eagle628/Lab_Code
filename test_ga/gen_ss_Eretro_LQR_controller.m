%% NOTE
% extend controller (\hat{K}*X Theorem 3.1)
% Rectifier is same at Simpel retrofit.
%%
classdef gen_ss_Eretro_LQR_controller < gen_ss
    
    properties
        model_local
        gen_ss
        params
        N
    end
    
    methods
        function obj = gen_ss_Eretro_LQR_controller(gen_ss, model_local)
            obj.gen_ss = gen_ss;
            obj.model_local = model_local;
            obj.params = gen_ss.params;
            obj.N = obj.gen_ss.N;
        end
        
        function [Anew, Bnew, Cnew, Dnew, dA, dB, dC, dD] = get_ss(obj, theta)
            if nargin < 2
                theta = obj.get_params();
            end
            [AE, BE, CE, DE, dAE, dBE, dCE, dDE] = obj.gen_ss.get_ss(theta);
            AL = obj.model_local.A;
            BL = obj.model_local(:, 'v').B;
            Cy = obj.model_local('y', :).C;
            Cw = obj.model_local('w', :).C;
            % calculate K
            [Anew, Bnew, Cnew] = make_new_local(obj.model_local, AE,BE,CE,DE);
            K = -dlqr(Anew, Bnew, eye(size(Anew, 1)), eye(size(Bnew, 2)));
            Ky = K(:, 1:size(AL, 1));
            Kw = K(:, size(AL, 1)+1:end);
            K = [Kw, Ky];
            nE = obj.gen_ss.n;
            nL = size(AL, 1);
            ny = size(Cy, 1);
            nw = size(Cw, 1);
            Anew = [AE, BE*Cw; BL*CE, AL+BL*DE*Cw];
            Bnew = [zeros(nE, ny), BE; zeros(nL, ny), BL*DE];
%             Cnew = [zeros(ny, nE), Cy; zeros(nw, nL), Cw];
%             Dnew = [eye(ny), zeros(ny, nw); zeros(nw, ny), eye(nw)];
            Cnew = [zeros(ny, nE), Cy; eye(nE), zeros(nE, nL)];
            Dnew = [eye(ny), zeros(ny, nw); zeros(nE, ny), zeros(nE, nw)];
            Cnew = K*Cnew;
            Dnew = K*Dnew;
            dA = cell(obj.N, 1);
            dB = cell(obj.N, 1);
            dC = cell(obj.N, 1);
            dD = cell(obj.N, 1);
            if nargout > 4
                for itr = 1 : numel(dAE)
                    dA{itr} = [dAE{itr}, dBE{itr}*Cw; BL*dCE{itr}, AL+BL*dDE{itr}*Cw];
                    dB{itr} = [zeros(nE, ny), dBE{itr}; zeros(nL, ny), BL*dDE{itr}];
%                     dC{itr} = K*[zeros(ny, nE), Cy; zeros(nw, nL), Cw];
%                     dD{itr} = K*[eye(ny), zeros(ny, nw); zeros(nw, ny), eye(nw)];
                    dC{itr} = K*[zeros(ny, nE), Cy; eye(nE), zeros(nE, nL)];
                    dD{itr} = K*[eye(ny), zeros(ny, nw); zeros(nE, ny), zeros(nE, nw)];
                end
            else
                return;
            end
        end
        
        function theta = get_params(obj)
            theta = obj.gen_ss.get_params();
            
        end
        
        function set_params(obj, theta)
            obj.gen_ss.set_params(theta);
        end
    end
end


%% local
function [Anew, Bnew, Cnew] = make_new_local(local, AE,BE,CE,DE)
    AL = local.A;
    Bu = local(:,'u').B;
    Bv = local(:,'v').B;
    Cy = local('y',:).C;
    Cw = local('w',:).C;
    Anew = [AE, BE*Cw; Bv*CE, Bv*DE*Cw+AL];
    Bnew = [tools.zeros(AE,Bu);Bu];
    Cnew = [tools.zeros(Cy,AE), Cy; tools.zeros(Cw,AE), Cw];
end
