classdef one_step_predictive_observer < State_Estimator
    % one step perdivtive observer :
    % Predict the next environmental internal state 
    % from the current observation and input using the model that has been pre-trained.
    % Args
    % A : environment model A
    % B : environment model B
    % C : environment model C
    % L : observer gain
    
    % NOTE : predict type
    properties
        A
        B
        C
        L
        model_A
        model_B
        model_C
    end
    
    methods
        function obj = one_step_predictive_observer(A, B, C, L)
            if nargin < 4
                Q = 1000*eye(size(A, 1));
                R = eye(size(C, 1));
                L = -dlqr(A', C', Q, R)';
            end
            obj.model_A = A;
            obj.model_B = B;
            obj.model_C = C;
            obj.A = A + L*C;
            obj.B = [B, -L];
            obj.C = eye(size(A));
            obj.L = L;
            obj.internal_state = zeros(size(A, 1), 1);
        end
        
        function belief_state = estimate(obj, varargin)
            narginchk(1, inf)
            belief_state = obj.internal_state; 
        end
        
        function update_internal_state(obj, y, u, varargin)
            narginchk(3, inf)
            obj.internal_state = obj.A*obj.internal_state + obj.B*[u;y]; 
        end
        
        function initialize(obj)
            obj.internal_state = zeros(size(obj.internal_state));
        end
        
        function [AL,BL,CL] = connect(obj, plant)
            [Ap,Bp,Cp,~]  = ssdata(plant);
            AL = [Ap, tools.zeros(Ap, obj.A);-obj.L*Cp, obj.A];
            BL = [Bp; obj.model_B];
            CL = [tools.zeros(obj.C, Ap), obj.C];
        end
    end
end
