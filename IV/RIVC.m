%% NOTE
% Identification of Continuous-time Models from Sampled Data.pdf 
% Capter 4 
% RIVC implement
% 
% Doubt :
% noise free output x is not using at Instrumental Variables.
% We need not to calculate "x".

% na,nb,nc,nd are not dimension but number of variables.
%% main function
function sys = RIVC(data,condition,lambda)
% % % % G = (s+1)*(s+5)/(s+3)/(s+10);
% % % % 
% % % % N = 10000;
% % % % Ts = 0.01;
% % % % 
% % % % G_d = c2d(G, Ts, 'foh');
% % % % H_d = tf([1,randn(1,kkk-1)],[1,randn(1,kkk-1)],Ts);
% % % % 
% % % % rng('shuffle')
% % % % u = randn(N,1)*1;
% % % % d = randn(N,1)*0.1;
% % % % y =  lsim(G_d, u) + lsim(H_d, d);
% % % % 
% % % % condition = [2,3,1,1];

    % load data
    y  = data.y;
    u  = data.u;
    Ts = data.Ts;
    N  = data.N;

    % coefiicinet number
    nf = condition(1); % nb (book)
    nb = condition(2); % na (book)
    nc = condition(3);
    nd = condition(4);
    % init system (prefiletring)
    s = tf('s');
    if nargin < 3
        lambda = 1;
    end
    init_sys = 1/(s+lambda)^nf;
    [num, den] = tfdata(init_sys,'v');
        
    
    for itr1 = 1 : 3
        x = lsim(c2d(tf(num, den), Ts, 'foh'), u);
        % Continuous-Time Solution eith Iterpolation by prefiltering
        y_fc = zeros(N, nf+1);
%         x_fc = zeros(N, nf+1);
        for itr2 = 0 : nf
            num = zeros(size(num));
            num(end-itr2) = 1;
            y_fc(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), y);
%             x_fc(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), x);
        end
        u_fc = zeros(N, nb);
        for itr2 = 0 : nb-1
            num = zeros(size(num));
            num(end-itr2) = 1;
            u_fc(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), u);
        end
        % noise system estimation
        if itr1 == 1
            noise_sys = tf(1,1,Ts);
            y_f = y_fc;
%             x_f = x_fc;
            u_f = u_fc;
        else
            xi = y-x;
            noise_sys = armax(iddata(xi,[],Ts),[nc,nd]);
            noise_sys = tf(noise_sys.A,noise_sys.C,Ts);
            % Discrete-Time Solution
            y_f = zeros(N, nf+1);
%             x_f = zeros(N, nf+1);
            for itr2 = 0 : nf
                y_f(:, end-itr2) = lsim(noise_sys, y_fc(:, end-itr2));
%                 x_f(:, end-itr2) = lsim(noise_sys, x_fc(:, end-itr2));
            end
            u_f = zeros(N, nb);
            for itr2 = 0 : nb-1
                u_f(:, end-itr2) = lsim(noise_sys, u_fc(:, end-itr2));
            end
        end
        % Instrumental variables
        phi = [-y_f(:,2:end), u_f];
        
        pre_rho = (phi'*phi)\(phi'*y_f(:,1));
        pre_rho = pre_rho';
        
        num = zeros(size(num));
        num(end-nb+1:end) = pre_rho(nf+1:nf+nb);
        den = [1,pre_rho(1:nf)];
        
%         sys = tf(num, den);
%         fbode(sys);
%         drawnow
    end
    sys = struct();
    sys.G = tf(num, den);
    sys.H = inv(noise_sys);
end


% bode(G,sys.G)
