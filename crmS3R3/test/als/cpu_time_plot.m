function cpu_time_plot(sols, varargin)

if nargin < 1
    error('Not enough arguments!')
else
    for j = 2:2:nargin
        switch lower(varargin{j-1})
            case {'error'}
                my_var = varargin{j};
            case{'combine'}
                combine = varargin{j};
            otherwise
                warning(['parameter ''' varargin{j-1} ''' not recognized']);
        end
    end
end
if ~exist('my_var','var')
    my_var = 'q';
end
if ~exist('combine','var')
    combine = 'no';
end

the_size = max(size(sols));
colors = linspecer(3);

for i = 1:the_size
    % separate figure per type of experiment
    switch combine
        case 'no'
            switch sols{i}.problem_name
                case 'roll-up'
                    if exist('n_roll_up','var')
                        figure(n_roll_up)
                        pl_rll.XData = [pl_rll.XData,sols{i}.err.abs.(my_var)];
                        pl_rll.YData = [pl_rll.YData,sols{i}.stats.cpu_time];
                    else
                        figure()
                        ax=gca;
                        ax.FontSize = 14;
                        set(gca, 'XScale', 'log', 'YScale', 'log')
                        n_roll_up = get(gcf,'Number');
                        pl_rll = loglog(sols{i}.err.abs.(my_var),sols{i}.stats.cpu_time,...
                            'Color',colors(1,:),'LineWidth',1.5);
                        title('Roll-up maneuver','FontSize',14)
                        xlabel(['Global error in ' my_var ],'FontSize',14)
                        ylabel('CPU time','FontSize',14)
                    end
                case 'flying_spaghetti'
                    if exist('n_flying_spaghetti','var')
                        figure(n_flying_spaghetti)
                        pl_fly.XData = [pl_fly.XData,sols{i}.err.abs.(my_var)];
                        pl_fly.YData = [pl_fly.YData,sols{i}.stats.cpu_time];
                    else
                        figure()
                        ax=gca;
                        ax.FontSize = 14;
                        set(gca, 'XScale', 'log', 'YScale', 'log')
                        n_flying_spaghetti = get(gcf,'Number');
                        pl_fly = loglog(sols{i}.err.abs.(my_var),sols{i}.stats.cpu_time,...
                            'Color',colors(2,:),'LineWidth',1.5);
                        title('Flying spaghetti','FontSize',14)
                        xlabel(['Global error in ' my_var ],'FontSize',14)
                        ylabel('CPU time','FontSize',14)
                    end
                otherwise
                    if exist('n_other','var')
                        figure(n_other)
                        pl_oth.XData = [pl_oth.XData,sols{i}.err.abs.(my_var)];
                        pl_oth.YData = [pl_oth.YData,sols{i}.stats.cpu_time];
                    else
                        figure()
                        ax=gca;
                        ax.FontSize = 14;
                        set(gca, 'XScale', 'log', 'YScale', 'log')
                        n_other = get(gcf,'Number');
                        pl_oth = loglog(sols{i}.err.abs.(my_var),sols{i}.stats.cpu_time,...
                            'Color',colors(3,:),'LineWidth',1.5);
                        title(strrep(sols{i}.problem_name,'_',' '),'FontSize',14)
                        xlabel(['Global error in ' my_var],'FontSize',14)
                        ylabel('CPU time','FontSize',14)
                    end
            end
            grid on
        case 'problem'
            if exist('n_figure','var')
                figure(n_figure)
            else
                figure()
                hold on
                n_figure = get(gcf,'Number');
                ax=gca;
                ax.FontSize = 14;
                set(gca, 'XScale', 'log', 'YScale', 'log')
                title('CPU time vs global error','FontSize',14)
                xlabel(['Global error in ' my_var ],'FontSize',14)
                ylabel('CPU time','FontSize',14)
                legend('Location','best','AutoUpdate','on')
                grid on
            end
            switch sols{i}.problem_name
                case 'roll-up'
                    if exist('pl_rll','var')
                        pl_rll.XData = [pl_rll.XData,sols{i}.err.abs.(my_var)];
                        pl_rll.YData = [pl_rll.YData,sols{i}.stats.cpu_time];
                    else
                        pl_rll = loglog(sols{i}.err.abs.(my_var),sols{i}.stats.cpu_time,...
                            'Color',colors(1,:),'LineWidth',1.5,'DisplayName','Roll-up');
                    end
                case 'flying_spaghetti'
                    if exist('pl_fly','var')
                        pl_fly.XData = [pl_fly.XData,sols{i}.err.abs.(my_var)];
                        pl_fly.YData = [pl_fly.YData,sols{i}.stats.cpu_time];
                    else
                        pl_fly = loglog(sols{i}.err.abs.(my_var),sols{i}.stats.cpu_time,...
                            'Color',colors(2,:),'LineWidth',1.5,'DisplayName','Flying spaghetti');
                    end
                otherwise
                    if exist('pl_oth','var')
                        pl_oth.XData = [pl_oth.XData,sols{i}.err.abs.(my_var)];
                        pl_oth.YData = [pl_oth.YData,sols{i}.stats.cpu_time];
                    else
                        pl_oth = loglog(sols{i}.err.abs.(my_var),sols{i}.stats.cpu_time,...
                            'Color',colors(3,:),'LineWidth',1.5,'DisplayName','Other test');
                    end
            end
        otherwise
            warning(['The case ''' combine ''' is not yet implemented!'])
    end

end

end