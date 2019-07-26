classdef test_apx_model < environment_model
    
    properties
        A
        B
        C
        D
        L
        Ts
        true_nx
        apx_nx
        nu
        ny
        true_sys
    end
    
    methods
        function obj = test_apx_model(Ts, true_dim, apx_dim)
            if nargin < 3
                apx_dim = 2;
            end
            if nargin < 2
                true_dim = 6;
            end
            if nargin < 1
                Ts = 0.01;
            end
            obj.true_nx = true_dim;
            obj.apx_nx = apx_dim;
            obj.Ts = Ts;
            sys = drss(true_dim);
            sys.Ts = Ts;
            obj.true_sys = sys;
            obj.nu = 1;
            obj.ny = 1;
            apx_sys = balred(sys, apx_dim);
            obj.A = apx_sys.A;
            obj.B = apx_sys.B;
            obj.C = apx_sys.C;
            obj.D = apx_sys.D;
        end
        
        function [ne_x, y] = dynamics(obj, pre_x, pre_u)
            ne_x = obj.true_sys.A*pre_x' + obj.true_sys.B*pre_u';
            y = obj.true_sys.C*pre_x' + obj.true_sys.D*pre_u';
        end
    end
end

