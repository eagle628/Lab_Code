classdef power_swing_2nd < handle
    % POWER_System Prameter Set
    %{
    m : Mass Parmeter
    d : Damper Parameter
    b : Input parameter
    name : object name 
    seed : Randam Number Seed
    %}
    properties
        m
        d
        b
        name
        seed
    end
    
    methods
        function obj = power_swing_2nd(m, d, b, name, seed)
            % Seed Control
            s = rng();
            if nargin == 5
                obj.seed = seed;
                rng(seed);
            else
                obj.seed = s;
            end
            % Nargin Missing
            if nargin < 4 || isempty(name)
                id = '';
                ss = rng();
                for itr = 1:6
                    rng('shuffle');
                    id = strcat(id, char(floor(double('a')+26*rand)));
                end
                rng(ss);
                name = id;
            end
            if nargin < 1 || isempty(m)
                m = 1+rand();
            end
            if nargin < 2 || isempty(d)
                d = (0.2+rand())*0.01;
            end
            if nargin < 3 || isempty(b)
                b = 1;
            end
            % Set Object Name
            obj.name = name;
            % Set Mass Parameter
            if numel(m) == 1
                obj.m = m;
            else
                obj.m = m(1)+rand()*(m(2)-m(1));
            end
            % Set Damper Parameter
            if numel(d) == 1
                obj.d = d;
            else
                obj.d = d(1)+rand()*(d(2)-d(1));
            end
            % Set Input
            obj.b = b;
        end
    end
end

