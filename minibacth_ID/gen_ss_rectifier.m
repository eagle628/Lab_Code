%% NOTE
% output rectifier form (2 input 1 output)
% input [w;v]
% output z (Filter y)
% asump : positive feedback (network)
%%
classdef gen_ss_rectifier < gen_ss
    
    properties
        model_local
        gen_ss
        params
        Cz
    end
    
    methods
        function obj = gen_ss_rectifier(gen_ss, model_local, Cz)
            obj.gen_ss = gen_ss;
            obj.model_local = model_local;
            obj.params = gen_ss.params;
            if nargin < 3
                Cz = eye(order(model_local));
            end
            obj.Cz = Cz;
        end
        
        function [Anew, Bnew, Cnew, Dnew, dA, dB, dC, dD] = get_ss(obj, theta)
            [Am, Bm, Cm, Dm, dAm, dBm, dCm, dDm] = obj.gen_ss.get_ss(theta);
            [A, L, G, ~] = ssdata(obj.model_local);
            Anew = [A+L*Dm*G, -L*Cm; -Bm*G, Am];% [local_state;env_state;]
            Bnew = [-L*Dm, L; Bm, tools.zeros(Bm, L)];
            %             Cnew = [-Dm*G, Cm];
            %             Dnew = [Dm, tools.zeros(Dm, L)];
            Cnew = [obj.Cz, tools.zeros(obj.Cz, Am)];
            Dnew = tools.zeros(obj.Cz, Bnew);
            dA = cell(numel(dAm), 1);
            dB = cell(numel(dAm), 1);
            dC = cell(numel(dAm), 1);
            dD = cell(numel(dAm), 1);
            for itr = 1:numel(dAm)
                dA{itr} = [L*dDm{itr}*G, -L*dCm{itr}; -dBm{itr}*G, dAm{itr}];
                dB{itr} = [-L*dDm{itr}, L*0; dBm{itr}, tools.zeros(Bm, L)];
                dC{itr} = tools.zeros(obj.Cz, Anew);
                dD{itr} = Dnew*0;
                %                 dC{itr} = [-dDm{itr}*G, dCm{itr}];
                %                 dD{itr} = [dDm{itr}, tools.zeros(Dm, L)];
            end
        end
        
        function theta = get_params(obj)
            theta = obj.gen_ss.get_params();
        end
        
        function set_params(obj, theta)
            obj.gen_ss.set_params(theta);
        end
        
        function varargout = con_blance(obj, theta)
           varargout = cell(nargout, 1);
           [varargout{:}] = obj.gen_ss.con_blance(theta);
        end
        
        function set_sys(obj, sys)
           obj.gen_ss.set_sys(sys); 
        end
        
    end
end

