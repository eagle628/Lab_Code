classdef figure_config
    
    
    properties
    end
    
    methods
    end
    
    methods(Static)
        %  Set Retro fiugre
        function H = set_figure_retro(name,position,title)
            if nargin < 1
                name = '';
            end
            if nargin < 2
                title =  {'y','\hat{y}','\tilde{y}'};
                position = [0,0,600,900];
            end
            if nargin < 3
                if iscell(position)
                    title = position;
                    position = [0,0,600,900];
                else
                    title =  {'y','\hat{y}','\tilde{y}'};
                end
            end
            H = struct();
            H.axes = struct();
            H.f = figure('Name',name,'position',position);
            H.axes.ax1 = subplot(3,1,1);
            H.axes.ax2 = subplot(3,1,2);
            H.axes.ax3 = subplot(3,1,3);
            H.axes.ax1.Title.Interpreter = 'latex';
            H.axes.ax1.Title.FontSize = 12;
            H.axes.ax1.Title.String = strcat('$',title{1},'$');
            H.axes.ax1.XGrid = 'on';
            H.axes.ax1.Box = 'on';
            H.axes.ax1.NextPlot = 'add';
            H.axes.ax2.Title.Interpreter = 'latex';
            H.axes.ax2.Title.FontSize = 12;
            H.axes.ax2.Title.String = strcat('$',title{2},'$');
            H.axes.ax2.XGrid = 'on';
            H.axes.ax2.Box = 'on';
            H.axes.ax2.NextPlot = 'add';
            H.axes.ax3.Title.FontSize = 12;
            H.axes.ax3.Title.Interpreter = 'latex';
            H.axes.ax3.Title.String = strcat('$',title{3},'$');
            H.axes.ax3.XGrid = 'on';
            H.axes.ax3.Box = 'on';
            H.axes.ax3.NextPlot = 'add';
        end
        
        % Plot Retro y & \hat{y}
        function ytilde = plot_retro1(axes,time,y,yhat,config)
            % Row nad Colums of y,ydesired is same.
            % times row is same(usually colums is 1.)
            % Premised set_figure_reto H.axes
            if nargin < 5
                config = {'linewidth',0.8};
            end
            % all y
            plot(axes.ax1,time,y,config{:});
            if ~isempty(yhat)
                % yhat
                plot(axes.ax2,time,yhat,config{:});
                % ytilde
                ytilde = y-yhat;
                plot(axes.ax3,time,ytilde,config{:});
            end
        end
        
        % Plot Retro y & \tilde{y}
        function yhat = plot_retro2(axes,time,y,ytilde,config)
            % Row nad Colums of y,ydesired is same.
            % times row is same(usually colums is 1.)
            % Premised set_figure_reto H.axes
            if nargin < 5
                config = {'linewidth',0.8};
            end
            % all y
            plot(axes.ax1,time,y,config{:});
            if ~isempty(ytilde)
                % yhat
                yhat = y-ytilde;
                plot(axes.ax2,time,yhat,config{:});
                % error xhat
                plot(axes.ax3,time,ytilde,config{:});
            end
        end
        % Set Bode Plot(SISO)
        function H = set_figure_bode(name,position)
            H = struct();
            H.axes = struct();
            if nargin < 2
                position = [0,0,600,600];
            end
            if nargin < 1
                name = '';
            end
            H.f = figure('Name',strcat('Bode_',name),'Position',position);
            H.axes.ax1 = subplot(2,1,1);
            hold on;
            H.axes.ax2 = subplot(2,1,2);
            hold on;

            H.axes.ax1.XScale = 'log';
            H.axes.ax1.XMinorGrid = 'on';
            H.axes.ax1.YMinorGrid = 'on';
            H.axes.ax1.Box = 'on';
            H.axes.ax1.XLabel.String = 'Frequency [rad/s]';
            H.axes.ax1.YLabel.String = 'Gain [dB]';

            H.axes.ax2.XScale = 'log';
            H.axes.ax2.XMinorGrid = 'on';
            H.axes.ax2.YMinorGrid = 'on';
            H.axes.ax2.Box = 'on';
            H.axes.ax2.XLabel.String = 'Frequency [rad/s]';
            H.axes.ax2.YLabel.String = 'Phase [deg]';
        end
        % Plot Bode
        function plot_bode(axes,omega,magnitude,phase,config)
            if nargin < 5
                config = {'linewidth',0.8};
            end
            semilogx(axes.ax1,omega,mag2db(magnitude),config{:});
            semilogx(axes.ax2,omega,phase,config{:});
        end
        
    end
    
end

