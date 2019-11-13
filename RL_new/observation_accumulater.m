classdef observation_accumulater < State_Estimator
    % Type that uses accumulated N observations as belief states.
    % Args
    % ny : RL environment observation number.
    % accumulate_N : Number to accumulate observation.
    % properties
    % A,B,C,D : LTI system parameter
    
    % NOTE : Filtered type
    properties
       A
       B
       C
       D
       ny
       accumulate_N
    end
    
    methods
        function obj = observation_accumulater(ny, accumulate_N)
            if nargin < 2
                accumulate_N = 1;
            end
            obj.ny = ny;
            obj.accumulate_N = accumulate_N;
            A = diag(ones(1,(accumulate_N-1)*ny),-ny);
            B = zeros(size(A,1),ny);
            B(1:ny,1:ny) = eye(ny);
            obj.A = A;
            obj.B = B;
            obj.C = A;
            obj.D = B;
            obj.internal_state = zeros(size(A, 1), 1);
        end
        
        function belief_state = estimate(obj, y)
            narginchk(2, inf)
            belief_state = obj.C*obj.internal_state + obj.D*y; 
        end
        
        function update_internal_state(obj, y, varargin)
            narginchk(2, inf)
            obj.internal_state = obj.A*obj.internal_state + obj.B*y; 
        end
        
        function initialize(obj)
            obj.internal_state = zeros(size(obj.internal_state));
        end
        
        function [AL,BL,CL] = connect(obj, plant)
            [Ap,Bp,Cp,~]  = ssdata(plant);
            AL = [Ap, tools.zeros(Ap, obj.A);obj.B*Cp, obj.A];
            BL = [Bp; tools.zeros(obj.A, Bp)];
            CL = [obj.D*Cp, obj.C];
        end
    end
end
