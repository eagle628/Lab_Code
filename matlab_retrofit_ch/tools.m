classdef tools
    
    
    properties
    end
    
    methods
    end
    
    methods(Static)
        function [w, ub, lb] = Lebesgue_sampler_equal(y, step)
            q = floor((y+step/2)./step);
            w = [1; diff(q)~=0];
            ub = (q+1)*step-step/2;
            lb = q*step-step/2;
        end
        
        function [w, ub, lb] = Lebesgue_sampler_thresholds(y, th)
            th = sort(th);
            if min(th) > min(y)
                th = [min(y)-abs(min(y))*0.1; th(:)];
            end
            if max(th) < max(y)
                th = [th(:); max(y)+abs(max(y))*0.1];
            end
            idx = 1:numel(th);
            q = floor(interp1(th, idx, y));
            w = [1; diff(q) ~= 0];
            lb = th(q);
            ub = th(q+1);
        end
        
        function [w, ub, lb] = Lebesgue_sampler(y, func)
            if isnumeric(func)
                if numel(func) == 1
                    [w, ub, lb] = tools.Lebesgue_sampler_equal(y, func);
                else
                    [w, ub, lb] = tools.Lebesgue_sampler_thresholds(y, func);
                end
            else
                
                
            end
            
        end
        
        function G = get_gaussian(x, mu, sigma2)
            G = exp(-(x-mu).^2/(2*sigma2))/sqrt(2*pi*sigma2);
        end
        
        function [ tfobj ] = sym2tf( symobj, vars, varvals)
            %SYM2TF Convert symbolic math rationals to transfer function
            % [ tfobj ] = sym2tf( symobj, vars, varvals)
            %     this function perform subs(symobj,vars,sym(varvals,'d'))
            
            if ~isa(symobj,'sym')
                tfobj=tf(symobj);
                return;
            end
            
            if nargin==3
                symobj=subs(symobj,vars,sym(varvals,'d'));
            end
            % if symobj==0
            %     tfobj=tf(0);
            %     return;
            % end
            if ~isa(symobj,'sym')
                tfobj=tf(symobj);
                return;
            end
            
            
            [n,d]=numden(symobj);
            dc=coeffs(d);
            n=n/dc(1);
            d=d/dc(1);
            num=sym2poly(n);
            den=sym2poly(d);
            fact=max(abs(den));
            num=num/fact;
            den=den/fact;
            
            tfobj=tf(num,den);
        end
        
        function [ tfobj ] = sym2dtf( symobj, vars, varvals, Ts)
            %SYM2TF Convert symbolic math rationals to transfer function
            % [ tfobj ] = sym2tf( symobj, vars, varvals)
            %     this function perform subs(symobj,vars,sym(varvals,'d'))
            
            if ~isa(symobj,'sym')
                tfobj=tf(symobj);
                return;
            end
            
            if nargin==4
                symobj=subs(symobj,vars,sym(varvals,'d'));
            end
            % if symobj==0
            %     tfobj=tf(0);
            %     return;
            % end
            if ~isa(symobj,'sym')
                tfobj=tf(symobj);
                return;
            end
            
            
            [n,d]=numden(symobj);
            dc=coeffs(d);
            n=n/dc(1);
            d=d/dc(1);
            num=sym2poly(n);
            den=sym2poly(d);
            fact=max(abs(den));
            num=num/fact;
            den=den/fact;
            
            tfobj=tf(num,den,Ts);
        end
        
        function [ p ] = pplot(x,y,z)
            %PPLOT Draw line as patch
            
            if nargin==2
                z=0;
            end
            
            x=reshape(x,[],1);
            y=reshape(y,[],1);
            
            p=patch([x;flipud(x)],[y;flipud(y)],z*ones(2*size(x,1),1),'b');
            
        end
        
        function close_plots(h)
            drawnow
            p1 = get(h(1), 'Position');
            p2 = get(h(2), 'Position');
            
            sp = 1-(p1(2) + p1(4)) + p2(2);
            set(h(1), 'Position',  p1 + [0,0,0,sp/4])
            set(h(2), 'Position', p2 + [0,0,0,sp/4+0.0936-0.1398/2])
        end
        
        function [varargout] = mybode2(rect, cell_sys, style)
            if nargin < 3 || isempty(style)
                style = 'b-';
            end
            if ~iscell(cell_sys)
                cell_sys = {cell_sys};
            end
            wmin = inf;
            wmax = -inf;
            for k=1:numel(cell_sys)
                [~,~,w] = bode(cell_sys{k});
                if wmin > min(w)
                    wmin = min(w);
                end
                if wmax < max(w)
                    wmax = max(w);
                end
            end
            nsys = numel(cell_sys);
            w = cell(nsys,1);
            MAG = cell(nsys,1);
            PHASE = cell(nsys,1);
            npoint = 0;
            for k=1:nsys
                [mag,phase,w{k}] = bode(cell_sys{k},{wmin, wmax});
                MAG{k} = squeeze(mag);
                phase = squeeze(phase);
                phase = phase - round(phase(1)/360)*360;
                PHASE{k} = phase;
                npoint = npoint + numel(w{k})+1;
            end
            if nargout < 2
                phase_p = nan(npoint, 1);
                mag_p = nan(npoint, 1);
                w_p = nan(npoint, 1);
                idx = 1;
                for k=1:nsys
                    phase_p(idx:idx+numel(PHASE{k})-1) = PHASE{k};
                    mag_p(idx:idx+numel(MAG{k})-1) = MAG{k};
                    w_p(idx:idx+numel(w{k})-1) = w{k};
                    idx = idx + numel(w{k}) + 1;
                end
                try
                    h = rect;
                    get(h(1), 'Position');
                    tools.myplot2(h(1), w_p, 20*log10(mag_p), style);
                    tools.myplot2(h(2), w_p, phase_p, style);
                catch
                    h = tools.mysubplot(w_p, {20*log10(mag_p), phase_p}, rect, style);
                    set(h(1), 'XScale', 'log')
                    set(h(2), 'XScale', 'log')
                    set(h(1),'xticklabel',[]);
                    xlabel(h(2), 'Frequency [rad/s]')
                    ylabel(h(1), 'Magnitude [dB]')
                    ylabel(h(2), 'Phase [deg]')
                    drawnow
                    p1 = get(h(1), 'Position');
                    p2 = get(h(2), 'Position');
                    
                    sp = 1-(p1(2) + p1(4)) + p2(2);
                    set(h(1), 'Position',  p1 + [0,0,0,sp/4])
                    set(h(2), 'Position', p2 + [0,0,0,sp/4+0.0936-0.1398/2])
                end
                varargout(1) = {h};
            else
                varargout(1)={MAG};
                varargout(2)={PHASE};
                varargout(3)={w};
            end
            
        end
        
        function [ varargout ] = mybode( varargin )
            fontName = 'Times New Roman';
            if nargin == 1&&iscell(varargin{1})
                varargin = varargin{1};
            end
            styles = {};
            if nargin == 2 && iscell(varargin{1}) && iscell(varargin{2})
                styles = varargin{2};
                varargin = varargin{1};
            end
            W = [];
            for k=1:numel(varargin)
                [~,~,w] = bode(varargin{k});
                W = [W;w];
            end
            w = sort(unique(W));
            MAG = zeros(numel(w),numel(varargin));
            PHASE = zeros(numel(w),numel(varargin));
            
            for k=1:numel(varargin)
                [mag,phase] = bode(varargin{k},w);
                MAG(:,k) = squeeze(mag);
                phase = squeeze(phase);
                phase = phase - round(phase(1)/360)*360;
                PHASE(:,k) = phase;
            end
            if nargout==0
                figure('Position',[100 100 760 570]);
                if isempty(styles)
                    h1 = subplot(2,1,1);semilogx(w,20*log10(MAG), 'b-');
                else
                    for itr = 1:size(MAG, 2)
                        if isstr(styles{itr})
                            h1 = subplot(2,1,1);semilogx(w,20*log10(MAG(:,itr)), styles{itr});
                        else
                            h1 = subplot(2,1,1);semilogx(w,20*log10(MAG(:,itr)), 'Color', styles{itr});
                        end
                        hold on
                    end
                    hold off
                end
                grid on
                %    xlim([w(1) w(end)]);
                ylabel('Magnitude [dB]','FontName',fontName,'FontSize',22);
                set(gca,'xticklabel',[]);
                set(gca,'Position',[0.15 0.584 0.775 0.314]);
                set(gca,'xtick',10.^(ceil(log10(w(1))):floor(log10(w(end)))));
                set( gca, 'FontName',fontName,'FontSize',20 );
                if isempty(styles)
                    h2 = subplot(2,1,2);semilogx(w,PHASE,'b-');
                else
                    for itr = 1:size(MAG, 2)
                        if isstr(styles{itr})
                            h2 = subplot(2,1,2);semilogx(w,PHASE(:,itr), styles{itr});
                        else
                            h2 = subplot(2,1,2);semilogx(w,PHASE(:,itr), 'Color', styles{itr});
                        end
                        hold on
                    end
                    hold off
                end
                grid on
                linkaxes([h1,h2],'x');
                set(gca,'Position',[0.15 0.19 0.775 0.314]);
                set(gca,'xtick',10.^(ceil(log10(w(1))):floor(log10(w(end)))));
                set( gca, 'FontName',fontName,'FontSize',20 );
                xlim([w(1) w(end)]);
                xlabel('Frequency [rad/s]','FontName',fontName,'FontSize',22);
                ylabel('Phase [deg]','FontName',fontName,'FontSize',22);
            end
            if nargout~=0
                varargout(1)={MAG};
                varargout(2)={PHASE};
                varargout(3)={w};
            end
        end
        
        function [v, dv, dv_num] = test_diff(func, theta, h, pararell)
            [v, dv] = func(theta);
            if nargin < 3
                h = 1e-8;
                pararell = false;
            end
            dv_num = zeros(numel(v), numel(theta));
            parfor_progress(numel(theta));
            if pararell
                parfor itr = 1:numel(theta)
                    theta1 = theta;
                    theta1(itr) = theta1(itr) + h;
                    v1 = func(theta1);
                    dv_num(:, itr) = (v1-v)/h;
                    parfor_progress();
                end
            else
                for itr = 1:numel(theta)
                    theta1 = theta;
                    theta1(itr) = theta1(itr) + h;
                    v1 = func(theta1);
                    dv_num(:, itr) = (v1-v)/h;
                    parfor_progress();
                end
            end
            parfor_progress(0);
            if numel(v) == 1
                dv_num = dv_num';
            end
        end
        
        function [v, dv, dv_num] = test_diff2(func, theta, h, pararell)
            [v, dv] = func(theta);
            if nargin < 3
                h = 1e-8;
                pararell = false;
            end
            parfor_progress(numel(theta));
            dv_num = [];
            if pararell
                parfor itr = 1:numel(theta)
                    theta1 = theta;
                    theta1(itr) = theta1(itr) + h;
                    v1 = func(theta1);
                    dv_num = [dv_num, (v1-v)/h];
                    parfor_progress();
                end
            else
                for itr = 1:numel(theta)
                    theta1 = theta;
                    theta1(itr) = theta1(itr) + h;
                    v1 = func(theta1);
                    dv_num = [dv_num, (v1-v)/h]; %#ok
                    parfor_progress();
                end
            end
            parfor_progress(0);
            if numel(v) == 1
                dv_num = dv_num';
            end
        end
        
        
        function X=ellipse(mu,S,q)
            if nargin < 3
                q = S;
                S = mu;
                mu = zeros(size(S,1), 1);
            end
            
            s2=diag(S);
            rho=S(1,2)/prod(sqrt(s2));
            
            D2=-2*log(1-q/100);
            p=1/(2*pi*sqrt(det(S)))*exp(-1/2*D2);
            C=-2*log(2*pi*sqrt(det(S))*p);
            t=(0:360).';
            P=1/sqrt(2)*[1 -1;1 1];
            tmpy1=cosd(t)*sqrt(C*(1+rho));
            tmpy2=sind(t)*sqrt(C*(1-rho));
            
            [R,~]=size(mu);
            if R==length(mu)
                tmp=mu;
            else
                tmp=mu.';
            end
            
            X=bsxfun(@plus,sqrt(diag(s2))*P*[tmpy1.';tmpy2.'],tmp).';
        end
        
        
        
        function rmse = RMSE(y1, y2)
            rmse = sqrt(mean((y1-y2).^2));
        end
        
        
        function df = num_diff(f)
            h = 1e-8;
            df =  @(x) (f(x + h) - f(x))./h;
        end
        
        function b = myisinteger(k)
            if abs(k-round(k)) < 10*eps
                b = true;
            else
                b = false;
            end
        end
        function [h, f1] = mysubplot_nodisp(x, y, rect, style, varargin)
            if nargin < 3
                rect = [500 450];
            end
            if nargin < 4
                style = '';
            end
            if ~iscell(style)
                style_org = style;
                style = cell(size(y));
                for itr = 1:numel(style)
                    style{itr} = style_org;
                end
            end
            fontName = 'Times New Roman';
            f1 = figure('Position', [0, 0, rect],  'visible', 'off');
            h = zeros(numel(y), 1);
            for itr = 1:numel(y)
                h(itr) = subplot(numel(y), 1, itr);
                %                 patchline(x, y{itr}, 'LineWidth', 1);
                plot(x, y{itr}, style{itr}, 'LineWidth', 1, varargin{:});
                set(gca, 'FontName',fontName,'FontSize',20 );
                grid on
                box on
            end
            linkaxes(h,'x');
        end
        function [h, f1] = mysubplot(x, y, rect, style, varargin)
            if nargin < 3 || isempty(rect)
                rect = [500 450];
            end
            if iscell(rect)
                xx = x;
                x = y;
                y = rect;
                rect = xx;
            end
            if nargin < 4
                style = '';
            end
            if ~iscell(style)
                style_org = style;
                style = cell(size(y));
                for itr = 1:numel(style)
                    style{itr} = style_org;
                end
            end
            fontName = 'Times New Roman';
            b = false;
            for itr = 1:numel(varargin)
                if ischar(varargin{itr})
                    if strcmpi('linewidth', varargin{itr})
                        b =true;
                    end
                end
            end
            if ~b
               varargin = [varargin, {'LineWidth', 1}]; 
            end
            f1 = figure('Position', [0, 0, rect]);
            h = zeros(numel(y), 1);
            for itr = 1:numel(y)
                h(itr) = subplot(numel(y), 1, itr);
                %                 patchline(x, y{itr}, 'LineWidth', 1);
                plot(x, y{itr}, style{itr}, varargin{:});
                set(gca, 'FontName',fontName,'FontSize',20 );
                grid on
                box on
            end
            linkaxes(h,'x');
        end
        
        function h = mysubplot_nolink(x, y, rect, style)
            if nargin < 3
                rect = [500 450];
            end
            if nargin < 4
                style = '';
            end
            if ~iscell(style)
                style_org = style;
                style = cell(size(y));
                for itr = 1:numel(style)
                    style{itr} = style_org;
                end
            end
            fontName = 'Times New Roman';
            f1 = figure('Position', [0, 0, rect]);
            h = zeros(numel(y), 1);
            for itr = 1:numel(y)
                h(itr) = subplot(numel(y), 1, itr);
                %                 patchline(x, y{itr}, 'LineWidth', 1);
                plot(x{itr}, y{itr}, style{itr}, 'LineWidth', 1);
                set(gca, 'FontName',fontName,'FontSize',20 );
                grid on
                box on
            end
        end
        
        function [f1, ax] = myplot(rect, varargin)
            if nargin == 1
                varargin = {rect};
                rect = [];
            end
            if isempty(rect)
                rect = [500 450];
            end
            fontName = 'Times New Roman';
            f1 = figure('Position', [0, 0, rect]);
            b = false;
            for itr = 1:numel(varargin)
                if ischar(varargin{itr})
                    if strcmpi('linewidth', varargin{itr})
                        b =true;
                    end
                end
            end
            if b
                plot(varargin{:})
            else
                plot(varargin{:}, 'linewidth', 1);
            end
            ax = gca;
            set(ax, 'FontName',fontName,'FontSize',20 );
            grid on
            box on
            
        end
        
        function f1 = myplot_nodisp(rect, varargin)
            if isempty(rect)
                rect = [500 450];
            end
            fontName = 'Times New Roman';
            f1 = figure('Position', [0, 0, rect], 'visible', 'off');
            plot(varargin{:}, 'linewidth', 1)
            set(gca, 'FontName',fontName,'FontSize',20 );
            grid on
            box on
            
        end
        
        function fix_plot(h)
            fontName = 'Times New Roman';
            set(h, 'FontName',fontName,'FontSize',20 );
            grid on
            box on
        end
        
        function h = myplot2(h, varargin)
            hold(h, 'on');
            fontName = 'Times New Roman';
            b = false;
            for itr = 1:numel(varargin)
                if ischar(varargin{itr})
                    if strcmpi('linewidth', varargin{itr})
                        b =true;
                    end
                end
            end
            if b
                plot(h, varargin{:});
            else
                plot(h, varargin{:}, 'linewidth', 1);
            end
            %             plot(h, varargin{:}, 'linewidth', 1)
            set(gca, 'FontName',fontName,'FontSize',20 );
            grid on
            box on
            
        end
        function l = myhline(h, y, last, varargin)
            if nargin < 3
                last = true;
            end
            hold(h, 'on');
            range = xlim(h);
            y = y(:)';
            yy = [y;y;nan(1,numel(y))];
            yy = yy(:);
            xx = repmat([range, nan], 1, numel(y))';
            b = false;
            for itr = 1:numel(varargin)
                if ischar(varargin{itr})
                    if strcmpi('linewidth', varargin{itr})
                        b =true;
                    end
                end
            end
            if b
                plot(h, xx, yy, varargin{:});
            else
                plot(h, xx, yy, varargin{:}, 'linewidth', 1);
            end
            hold(h, 'off');
            g = get(h);
            l = g.Children(1);
            if last
                set(h, 'Children', [g.Children(2:end); g.Children(1)]);
            end
        end
        function l = myvline(h, x, last, varargin)
            if nargin < 3
                last = true;
            end
            hold(h, 'on');
            range = ylim(h);
            y = x(:)';
            yy = [y;y;nan(1,numel(y))];
            yy = yy(:);
            xx = repmat([range, nan], 1, numel(y))';
            plot(h, yy, xx, varargin{:});
            hold(h, 'off');
            g = get(h);
            l = g.Children(1);
            ylim(h, range);
            if last
                set(h, 'Children', [g.Children(2:end); g.Children(1)]);
                
            end
        end
        
        function sys = myrss(deg, Ts)
            sys = rss(deg);
            y = lsim(sys, ones(2,1), [0; Ts]);
            [msgstr, msgid] = lastwarn('');
            sys_min = minreal(sys);
            while ~strcmp(msgid, '') || size(sys_min.a,1) ~= deg || sum(pole(sys)==0)~=0
                sys = rss(deg);
                y = lsim(sys, ones(2,1), [0; Ts]);
                [msgstr, msgid] = lastwarn('');
                sys_min = minreal(sys);
            end
        end
        
        function start_cluster(n)
            c = parcluster('local');
            c.NumWorkers = n;
            parpool(c, c.NumWorkers);
        end
        
        function v = select_value(isA, a, b)
            if isA
                v = a;
            else
                v = b;
            end
        end
        
        function i = argmin(A)
            [~, i] = min(A);
        end
        function i = argmax(A)
            [~, i] = max(A);
        end
        function z = zeros(A, B)
            z = zeros(size(A, 1), size(B, 2));
        end
    end
    
end

