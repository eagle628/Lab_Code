%% NOTE
% Identification of Continuous-time Models from Sampled Data.pdf 
% Chapter 5
% section 5.5.3
% Two Step CLIVC2 algorithm
%% main
function G2 = CLIVC2(data,C,condition,lambda)
    % load data
    u  = data.u;
    y  = data.y;
    Ts = data.Ts;
    N  = length(u);

    % conditon
    na = condition(1);
    nb = condition(2);
    
    % contoroller
    if C.Ts == 0
        C_d = c2d(C, Ts, 'foh');
    else
        C_d = C;
    end
    % clculate r
    r = u + lsim(C_d, y);
    % chose CTfilter : f_c(p)
    if nargin < 4
        lambda = 1;
    end
    beta = 1;
    s = tf('s');
    init_sys = (beta/(s+lambda))^na;
    [num, den] = tfdata(init_sys,'v');

    for itr1 = 1:1

        % step : 1
        y_fc = zeros(N, na+1);
        for itr2 = 0 : na
            num = zeros(size(num));
            num(end-itr2) = 1;
            y_fc(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), y);
        end
        u_fc = zeros(N, nb);
        for itr2 = 0 : nb-1
            num = zeros(size(num));
            num(end-itr2) = 1;
            u_fc(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), u);
        end
        %  eq 5.38
        phi_fc = [-y_fc(:,2:end), u_fc];
        rho_1 = (phi_fc'*phi_fc)\(phi_fc'*y_fc(:,1));
        % step : 2
        G1 = tf((rho_1(na+1:end))',[1,(rho_1(1:na))']);
        G1_d = c2d(G1, Ts, 'foh');
        cloop_d_1 = loopsens(G1_d,C_d);
%         [num, den] = tfdata(G1,'v');
% % % %         if flag
% % % %             % Book Method
% % % %             % prefiltering r
% % % %             r_fc = zeros(N, max(na,nb));
% % % %             for itr2 = 0 : max(na,nb)-1
% % % %                 num = zeros(size(num));
% % % %                 num(end-itr2) = 1;
% % % %                 r_fc(:,end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), r);
% % % %             end
% % % %             % Eq 5.39
% % % %             x_fc = zeros(N, na);
% % % %             for itr2 = 0 : na-1
% % % %                 x_fc(:, end-itr2) = lsim(G1_d*cloop_d_1.Si, r_fc(:, end-itr2));
% % % %             end
% % % %             % Eq 5.40
% % % %             nu_fc = zeros(N, nb);
% % % %             for itr2 = 0 : nb-1
% % % %                 nu_fc(:, end-itr2) = lsim(cloop_d_1.Si, r_fc(:,end-itr2));
% % % %             end
% % % %         else
            % calculate noise free
            x = lsim(G1_d*cloop_d_1.Si, r);
            nu = lsim(cloop_d_1.Si, r);
            % Eq 5.39
            x_fc = zeros(N, na);
            for itr2 = 0 : na-1
                num = zeros(size(num));
                num(end-itr2) = 1;
                x_fc(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), x);
            end
            % Eq 5.40
            nu_fc = zeros(N, nb);
            for itr2 = 0 : nb-1
                num = zeros(size(num));
                num(end-itr2) = 1;
                nu_fc(:, end-itr2) = lsim(c2d(tf(num, den), Ts, 'foh'), nu);
            end
% % % %         end
        % Eq 5.41
        zeta_fc = [-x_fc, nu_fc];
        % Eq 5.42
        rho2 = (zeta_fc'*phi_fc)\(zeta_fc'*y_fc(:,1));
        G2 = tf((rho2(na+1:end))',[1,(rho2(1:na))']);
%         [num, den] = tfdata(G2,'v');

    end
end
