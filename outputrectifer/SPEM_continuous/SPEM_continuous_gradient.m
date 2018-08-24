%% Note
% Stabilized Prediciton Error Method
% Calculate continuous Gradient 
%% Main
function dPhat = SPEM_continuous_gradient(dim,K)
    K = tf2sym(K);

    a = sym('a',[dim, 1]);
    b = sym('b',[dim+1, 1]);
    
    syms s

    Phat = s.^(dim:-1:0)*b/(s.^(dim:-1:0)*[1;a]);
    Shat = 1/(1+Phat.*K);
    dPhat = cell( size(a,1)+size(b,1), 2);
    % input y
    for itr1 = 1 : size(a,1)
        dPhat{itr1,1} = diff(Shat,a(itr1));
    end
    for itr2 = 1 : size(b,1)
        dPhat{itr1+itr2,1} = diff(Shat,b(itr2));
    end
    % input u
    G = -Shat.*Phat;
    for itr1 = 1 : size(a,1)
        dPhat{itr1,2} = diff(G,a(itr1));
    end
    for itr2 = 1 : size(b,1)
        dPhat{itr1+itr2,2} = diff(G,b(itr2));
    end
end

%% local
function sym_sys = tf2sym(tf_sys)
    [ num ,den] = tfdata(tf_sys, 'v');
    dim = numel(num)-1; % Order
    d = sym('d',[numel(den), 1]); % Denominator parameter number
    n = sym('n',[numel(num), 1]); % Numerator parameter number
    
    syms s
    G = s.^(dim:-1:0)*n/(s.^(dim:-1:0)*d);
    sym_sys = subs(G, [d ; n], [den, num]');
end
    