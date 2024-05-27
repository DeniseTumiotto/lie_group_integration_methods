function orders_plot(sols, varargin)

    % if nargin < 1
    %     error('Not enough arguments!')
    % else
    %     for j = 2:2:nargin
    %         switch lower(varargin{j-1})
    %             case {'variable'}
    %                 my_var = varargin{j};
    %             case {'side'}
    %                 side = varargin{j};
    %             otherwise
    %                 warning(['parameter ''' varargin{j-1} ''' not recognized']);
    %         end
    %     end
    % end
    % if ~exist('my_var','var')
    %     my_var = 'all';
    % end
    % if ~exist('side','var')
    %     side = 0;
    % end
    
    % As for start, I will plot err.abs.q and err.abs.l side by side!
    the_size = max(size(sols));
    colors = linspecer(5);
    spaces = [];
    orders = [];

    for i = 1:the_size
        if sols{i}.err.abs.q == 0 || any(isnan([sols{i}.err.abs.q, sols{i}.err.abs.v, sols{i}.err.abs.l]))
            continue
        end
        j = find(spaces==sols{i}.problemset);
        if isempty(j)
            spaces = [spaces sols{i}.problemset];
            figure(length(spaces))
            sgtitle(num2str(sols{i}.problemset))
            k = find(orders==sols{i}.order);
            if isempty(k)
                orders = [orders sols{i}.order];
                subplot(1,2,1)
                hold on
                ax=gca;
                ax.FontSize = 14;
                ax.PlotBoxAspectRatio = [1 1 1];
                ax.XLim = [1e-4 1e-2];
                ax.YLim = [1e-10 1e0];
                ax.XTick = [1e-4 1e-3 1e-2];
                ax.XTickLabel = {'10^{-4}' '10^{-3}' '10^{-2}'};
                ax.YTick = [1e-10 1e-8 1e-6 1e-4 1e-2 1e0];
                ax.YTickLabel = {'10^{-10}', '10^{-8}', '10^{-6}', '10^{-4}', '10^{-2}', '10^{0}'};
                set(gca, 'XScale', 'log', 'YScale', 'log')
                % axis([10^-4 10^-2 10^-15 10^-4])
                title('\bf{Convergence order: $q$}','FontSize',14, 'Interpreter', 'latex')
                xlabel('Time step size (h)','FontSize',14, 'Interpreter', 'latex')
                ylabel('Global error','FontSize',14, 'Interpreter', 'latex')
                legend('Location','southeast','AutoUpdate','on')
                grid on
                plt_q{length(spaces),length(orders)} = loglog(sols{i}.h, sols{i}.err.abs.q, ...
                'Color',colors( - 1 + sols{i}.order,:),'LineWidth',1.5, ...
                'DisplayName',['p = ' num2str(sols{i}.order)]);
                subplot(1,2,2)
                hold on
                ax=gca;
                ax.FontSize = 14;
                ax.PlotBoxAspectRatio = [1 1 1];
                ax.XLim = [1e-4 1e-2];
                ax.YLim = [1e-8 1e2];
                ax.XTick = [1e-4 1e-3 1e-2];
                ax.XTickLabel = {'10^{-4}', '10^{-3}', '10^{-2}'};
                ax.YTick = [1e-8 1e-6 1e-4 1e-2 1e0 1e2];
                ax.YTickLabel = {'10^{-8}', '10^{-6}', '10^{-4}', '10^{-2}', '10^{0}', '10^{2}'};
                set(gca, 'XScale', 'log', 'YScale', 'log')
                % axis([10^-4 10^-2 10^-15 10^-4])
                title('\bf{Convergence order: $\lambda$}','FontSize',14, 'Interpreter', 'latex')
                xlabel('Time step size (h)','FontSize',14, 'Interpreter', 'latex')
                ylabel('Global error','FontSize',14, 'Interpreter', 'latex')
                legend('Location','southeast','AutoUpdate','on')
                grid on
                plt_l{length(spaces),length(orders)} = loglog(sols{i}.h, sols{i}.err.abs.l, ...
                'Color',colors( - 1 + sols{i}.order,:),'LineWidth',1.5, ...
                'DisplayName',['p = ' num2str(sols{i}.order)]);
            else
                subplot(1,2,1)
                hold on
                ax=gca;
                ax.FontSize = 14;
                ax.PlotBoxAspectRatio = [1 1 1];
                ax.XLim = [1e-4 1e-2];
                ax.YLim = [1e-10 1e0];
                ax.XTick = [1e-4 1e-3 1e-2];
                ax.XTickLabel = {'10^{-4}', '10^{-3}', '10^{-2}'};
                ax.YTick = [1e-10 1e-8 1e-6 1e-4 1e-2 1e0];
                ax.YTickLabel = {'10^{-10}', '10^{-8}', '10^{-6}', '10^{-4}', '10^{-2}', '10^{0}'};
                set(gca, 'XScale', 'log', 'YScale', 'log')
                % axis([10^-4 10^-2 10^-15 10^-4])
                title('\bf{Convergence order: $q$}','FontSize',14, 'Interpreter', 'latex')
                xlabel('Time step size (h)','FontSize',14, 'Interpreter', 'latex')
                ylabel('Global error','FontSize',14, 'Interpreter', 'latex')
                legend('Location','southeast','AutoUpdate','on')
                grid on
                plt_q{length(spaces),k} = loglog(sols{i}.h, sols{i}.err.abs.q, ...
                'Color',colors( - 1 + sols{i}.order,:),'LineWidth',1.5, ...
                'DisplayName',['p = ' num2str(sols{i}.order)]);
                subplot(1,2,2)
                hold on
                ax=gca;
                ax.FontSize = 14;
                ax.PlotBoxAspectRatio = [1 1 1];
                ax.XLim = [1e-4 1e-2];
                ax.YLim = [1e-8 1e2];
                ax.XTick = [1e-4 1e-3 1e-2];
                ax.XTickLabel = {'10^{-4}', '10^{-3}', '10^{-2}'};
                ax.YTick = [1e-8 1e-6 1e-4 1e-2 1e0 1e2];
                ax.YTickLabel = {'10^{-8}', '10^{-6}', '10^{-4}', '10^{-2}', '10^{0}', '10^{2}'};
                set(gca, 'XScale', 'log', 'YScale', 'log')
                % axis([10^-4 10^-2 10^-15 10^-4])
                title('\bf{Convergence order: $\lambda$}','FontSize',14, 'Interpreter', 'latex')
                xlabel('Time step size (h)','FontSize',14, 'Interpreter', 'latex')
                ylabel('Global error','FontSize',14, 'Interpreter', 'latex')
                legend('Location','southeast','AutoUpdate','on')
                grid on
                plt_l{length(spaces),k} = loglog(sols{i}.h, sols{i}.err.abs.l, ...
                'Color',colors( - 1 + sols{i}.order,:),'LineWidth',1.5, ...
                'DisplayName',['p = ' num2str(sols{i}.order)]);
            end
        else
            figure(j)
            k = find(orders==sols{i}.order);
            if isempty(k)
                orders = [orders sols{i}.order];
                subplot(1,2,1)
                hold on
                ax=gca;
                ax.FontSize = 14;
                ax.PlotBoxAspectRatio = [1 1 1];
                ax.XLim = [1e-4 1e-2];
                ax.YLim = [1e-10 1e0];
                ax.XTick = [1e-4 1e-3 1e-2];
                ax.XTickLabel = {'10^{-4}', '10^{-3}', '10^{-2}'};
                ax.YTick = [1e-10 1e-8 1e-6 1e-4 1e-2 1e0];
                ax.YTickLabel = {'10^{-10}', '10^{-8}', '10^{-6}', '10^{-4}', '10^{-2}', '10^{0}'};
                set(gca, 'XScale', 'log', 'YScale', 'log')
                % axis([10^-4 10^-2 10^-15 10^-4])
                title('\bf{Convergence order: $q$}','FontSize',14, 'Interpreter', 'latex')
                xlabel('Time step size (h)','FontSize',14, 'Interpreter', 'latex')
                ylabel('Global error','FontSize',14, 'Interpreter', 'latex')
                legend('Location','southeast','AutoUpdate','on')
                grid on
                plt_q{j,length(orders)} = loglog(sols{i}.h, sols{i}.err.abs.q, ...
                'Color',colors( - 1 + sols{i}.order,:),'LineWidth',1.5, ...
                'DisplayName',['p = ' num2str(sols{i}.order)]);
                subplot(1,2,2)
                hold on
                ax=gca;
                ax.FontSize = 14;
                ax.PlotBoxAspectRatio = [1 1 1];
                ax.XLim = [1e-4 1e-2];
                ax.YLim = [1e-8 1e2];
                ax.XTick = [1e-4 1e-3 1e-2];
                ax.XTickLabel = {'10^{-4}', '10^{-3}', '10^{-2}'};
                ax.YTick = [1e-8 1e-6 1e-4 1e-2 1e0 1e2];
                ax.YTickLabel = {'10^{-8}', '10^{-6}', '10^{-4}', '10^{-2}', '10^{0}', '10^{2}'};
                set(gca, 'XScale', 'log', 'YScale', 'log')
                % axis([10^-4 10^-2 10^-15 10^-4])
                title('\bf{Convergence order: $\lambda$}','FontSize',14, 'Interpreter', 'latex')
                xlabel('Time step size (h)','FontSize',14, 'Interpreter', 'latex')
                ylabel('Global error','FontSize',14, 'Interpreter', 'latex')
                legend('Location','southeast','AutoUpdate','on')
                grid on
                plt_l{j,length(orders)} = loglog(sols{i}.h, sols{i}.err.abs.l, ...
                'Color',colors( - 1 + sols{i}.order,:),'LineWidth',1.5, ...
                'DisplayName',['p = ' num2str(sols{i}.order)]);
            else
                if isempty(plt_q{j,k})
                    subplot(1,2,1)
                    hold on
                    ax=gca;
                    ax.FontSize = 14;
                    ax.PlotBoxAspectRatio = [1 1 1];
                    ax.XLim = [1e-4 1e-2];
                    ax.YLim = [1e-10 1e0];
                    ax.XTick = [1e-4 1e-3 1e-2];
                    ax.XTickLabel = {'10^{-4}', '10^{-3}', '10^{-2}'};
                    ax.YTick = [1e-10 1e-8 1e-6 1e-4 1e-2 1e0];
                    ax.YTickLabel = {'10^{-10}', '10^{-8}', '10^{-6}', '10^{-4}', '10^{-2}', '10^{0}'};
                    set(gca, 'XScale', 'log', 'YScale', 'log')
                    % axis([10^-4 10^-2 10^-15 10^-4])
                    title('\bf{Convergence order: $q$}','FontSize',14, 'Interpreter', 'latex')
                    xlabel('Time step size (h)','FontSize',14, 'Interpreter', 'latex')
                    ylabel('Global error','FontSize',14, 'Interpreter', 'latex')
                    legend('Location','southeast','AutoUpdate','on')
                    grid on
                    plt_q{j,k} = loglog(sols{i}.h, sols{i}.err.abs.q, ...
                    'Color',colors( - 1 + sols{i}.order,:),'LineWidth',1.5, ...
                    'DisplayName',['p = ' num2str(sols{i}.order)]);
                    subplot(1,2,2)
                    hold on
                    ax=gca;
                    ax.FontSize = 14;
                    ax.PlotBoxAspectRatio = [1 1 1];
                    ax.XLim = [1e-4 1e-2];
                    ax.YLim = [1e-8 1e2];
                    ax.XTick = [1e-4 1e-3 1e-2];
                    ax.XTickLabel = {'10^{-4}', '10^{-3}', '10^{-2}'};
                    ax.YTick = [1e-8 1e-6 1e-4 1e-2 1e0 1e2];
                    ax.YTickLabel = {'10^{-8}', '10^{-6}', '10^{-4}', '10^{-2}', '10^{0}', '10^{2}'};
                    set(gca, 'XScale', 'log', 'YScale', 'log')
                    % axis([10^-4 10^-2 10^-15 10^-4])
                    title('\bf{Convergence order: $\lambda$}','FontSize',14, 'Interpreter', 'latex')
                    xlabel('Time step size (h)','FontSize',14, 'Interpreter', 'latex')
                    ylabel('Global error','FontSize',14, 'Interpreter', 'latex')
                    legend('Location','southeast','AutoUpdate','on')
                    grid on
                    plt_l{j,k} = loglog(sols{i}.h, sols{i}.err.abs.l, ...
                    'Color',colors( - 1 + sols{i}.order,:),'LineWidth',1.5, ...
                    'DisplayName',['p = ' num2str(sols{i}.order)]);
                else
                    plt_q{j,k}.XData = [plt_q{j,k}.XData sols{i}.h];
                    plt_q{j,k}.YData = [plt_q{j,k}.YData sols{i}.err.abs.q];
                    plt_l{j,k}.XData = [plt_l{j,k}.XData sols{i}.h];
                    plt_l{j,k}.YData = [plt_l{j,k}.YData sols{i}.err.abs.l];
                end
            end
        end

    end

end