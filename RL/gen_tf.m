classdef gen_tf < handle
    %UNTITLED3 このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        params
        theta
        n % tf order
        m % input number
        l % output number
        N % parmeter length
    end
    
    methods
        function obj = gen_tf(n, m, l)
            obj.n = n;
            obj.m = m;
            obj.l = l;
            obj.N = n+(n+1)*(m*l);
            obj.params = cell(obj.N, 1);           
            obj.theta = zeros(obj.N, 1);
            idx = 1;
            for itr1 = n-1:-1:0
                obj.params{idx} = sprintf('theta_den_all_%d', itr1);
                idx = idx + 1;
            end
            for itr1 = 1 : l
                for itr2 =  1 : m
                    for itr3 = n:-1:0
                        obj.params{idx} = sprintf('theta_num_h%d%d_%d', itr1, itr2, itr3);
                        idx = idx + 1;
                    end
                end
            end
            
        end
        
        function set_params(obj, theta)
            obj.theta = theta;
        end
        
        function theta = get_params(obj)
            theta = obj.theta;
        end
        
        function set_sys(obj, sys)
            [num, den] = tfdata(sys);
            den = den{1,1};
            obj.theta(1:obj.n, :) = den(2:end)';
            idx = obj.n+1;
            for itr1 = 1 : obj.l
                for itr2 =  1 : obj.m
                    obj.theta(idx:idx+(obj.n), 1) = num{itr1, itr2}';
                    idx = idx + (obj.n) + 1;
                end
            end
        end
        
        function [A,B,C,D,dA,dB,dC,dD] = get_sys(obj, theta)
            if nargin < 2
                theta = obj.get_params();
            end
            den = {[1, theta(1:obj.n, 1)']};
            den_set = repmat(den, obj.l, obj.m);
            num_set = cell(obj.l, obj.m);
            idx1 = obj.n+1;
            for itr1 = 1 : obj.l
                for itr2 = 1 : obj.m
                    num_set{itr1, itr2} = theta(idx1:idx1+obj.n, :)';
                    idx1 = idx1 + (obj.n) + 1;
                end
            end
            sys = tf(num_set, den_set);
            if nargout == 1
                A = sys;
                return
            end
            [A,B,C,D] = ssdata(sys);
            if nargout > 2
                dA = cell(obj.N, 1);
                dB = cell(obj.N, 1);
                dC = cell(obj.N, 1);
                dD = cell(obj.N, 1);
                idx1 = 1;
                % dennominator
                % 商の微分の一部のみ（元のシステムと直列に接続する必要あり．）
                for itr1 = 1 : obj.n
                    num = zeros(1, 1+obj.n);
                    num(itr1) = -1;
                    sys = tf(repmat({num},1,obj.l), den_set(:,1)');
                    [dA{idx1,1},dB{idx1,1},dC{idx1,1},dD{idx1,1}] = ssdata(sys);
                    idx1 = idx1 + 1;
                end
                % numerator
                for itr1 = 1 : obj.l
                    for itr2 = 1 : obj.m
                        for itr3 = 1 : obj.n+1
                             num = zeros(1, 1+obj.n);
                             num_set = repmat({num}, obj.l, obj.m);
                             num_set{itr1, itr2}(itr3) = 1;
                             sys = tf(num_set, den_set);
                             [dA{idx1,1},dB{idx1,1},dC{idx1,1},dD{idx1,1}] = ssdata(sys);
                             idx1 = idx1 + 1;
                        end
                    end
                end
            end
        end
    end
end

