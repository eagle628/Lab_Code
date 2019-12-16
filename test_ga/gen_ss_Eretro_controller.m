%% NOTE
% extend controller (\hat{K}*X Theorem 3.1)
% Rectifier is same at Simpel retrofit.
%%
classdef gen_ss_Eretro_controller < gen_ss
    
    properties
        model_local
        gen_ss
        params
        K
        N
        n
    end
    
    methods
        function obj = gen_ss_Eretro_controller(gen_ss, model_local, K)
            obj.gen_ss = gen_ss;
            obj.model_local = model_local;
            obj.params = gen_ss.params;
            if nargin < 3
                K = randn(...
                    length(model_local.InputGroup.u),...
                    length(model_local.OutputGroup.y)+length(model_local.OutputGroup.w) ...
                    );
            end
            obj.K = K;
            obj.N = obj.gen_ss.N+numel(obj.K);
            obj.n = obj.gen_ss.n+order(model_local);
        end
        
        function [Anew, Bnew, Cnew, Dnew, dA, dB, dC, dD] = get_ss(obj, theta)
            if nargin < 2
                theta = obj.get_params();
            end
            ss_theta = theta(1:end-numel(obj.K));
            K_theta = reshape(theta(end-numel(obj.K)+1:end), size(obj.K,1), []);
            [AE, BE, CE, DE, dAE, dBE, dCE, dDE] = obj.gen_ss.get_ss(ss_theta);
            AL = obj.model_local.A;
            BL = obj.model_local(:, 'v').B;
            Cy = obj.model_local('y', :).C;
            Cw = obj.model_local('w', :).C;
            
            nE = obj.gen_ss.n;
            nL = size(AL, 1);
            ny = size(Cy, 1);
            nw = size(Cw, 1);
            Anew = [AE, BE*Cw; BL*CE, AL+BL*DE*Cw];
            Bnew = [zeros(nE, ny), BE; zeros(nL, ny), BL*DE];
            Cnew = [zeros(ny, nE), Cy; zeros(nw, nL), Cw];
            Dnew = [eye(ny), zeros(ny, nw); zeros(nw, ny), eye(nw)];
            Cnew = K_theta*Cnew;
            Dnew = K_theta*Dnew;
            dA = cell(obj.N, 1);
            dB = cell(obj.N, 1);
            dC = cell(obj.N, 1);
            dD = cell(obj.N, 1);
            if nargout > 4
                for itr = 1 : numel(dAE)
                    dA{itr} = [dAE{itr}, dBE{itr}*Cw; BL*dCE{itr}, AL+BL*dDE{itr}*Cw];
                    dB{itr} = [zeros(nE, ny), dBE{itr}; zeros(nL, ny), BL*dDE{itr}];
                    dC{itr} = K_theta*[zeros(ny, nE), Cy; zeros(nw, nL), Cw];
                    dD{itr} = K_theta*[eye(ny), zeros(ny, nw); zeros(nw, ny), eye(nw)];
                end
                b = itr;
                AE_ = AE*0;
                BE_ = BE*0;
                CE_ = CE*0;
                DE_ = DE*0;
                for itr = 1 : numel(obj.K)
                    dA{b+itr} = [AE_, BE_*Cw; BL*CE_, AL+BL*DE_*Cw];
                    dB{b+itr} = [zeros(nE, ny), BE_; zeros(nL, ny), BL*DE_];
                    dK = zeros(1, numel(obj.K));
                    dK(itr) = 1;
                    dK = reshape(dK, size(obj.N, 1), []);
                    dC{b+itr} = dK*[zeros(ny, nE), Cy; zeros(nw, nL), Cw];
                    dD{b+itr} = dK*[eye(ny), zeros(ny, nw); zeros(nw, ny), eye(nw)];
                end
            else
                return;
            end
        end
        
        function theta = get_params(obj)
            theta = [obj.gen_ss.get_params(), reshape(obj.K, 1, [])];            
        end
        
        function set_params(obj, theta)
            Ksize = numel(obj.K);
            obj.gen_ss.set_params(theta(1:end-Ksize));
            obj.K = reshape(theta(end-Ksize+1:end), size(obj.K, 1), []);
        end
    end
end

