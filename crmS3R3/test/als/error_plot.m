function error_plot(sols)

    the_size = max(size(sols));
    colors = linspecer(5);
    the_name = 'n';t_end = 1e7;
    figure()
    
    for i = 1:the_size
        n_t = max(size(sols{i}.rslt.t));
        for j = 1:n_t
            if sols{i}.rslt.t(j) > t_end
                break
            end
        end
        time_interval = sols{i}.rslt.t(2:j);
        my_problem = strrep(string(sols{i}.(the_name)),'_',' ');
        if size(my_problem,2)>1
            my_problem = my_problem(1);
        end
        if strcmp(sols{i}.problem_name,'roll-up')
            subplot(1,2,1)
            ax=gca;
            ax.FontSize = 14;
            ax.PlotBoxAspectRatio = [1 1 1];
            hold on
            title('Roll-up maneuver','FontSize',14)
            xlabel('Time (t)',...
            'FontSize',14,'Interpreter','latex')
            ylabel('$h_{\textrm{new}}$',...
            'FontSize',14,'Interpreter','latex')
            grid on
        else
            subplot(1,2,2)
            ax=gca;
            ax.FontSize = 14;
            ax.PlotBoxAspectRatio = [1 1 1];
            hold on
            title('Flying spaghetti','FontSize',14)
            xlabel('Time (t)',...
            'FontSize',14,'Interpreter','latex')
            ylabel('$h_{\textrm{new}}$',...
            'FontSize',14,'Interpreter','latex')
            grid on
        end
        if sols{i}.n == 8
            col = 3;
        else
            col = 4;
        end
        plot(time_interval,sols{i}.h(1:j-1),...
        'LineWidth',1.5,'DisplayName',['N = ' num2str(my_problem) ', CPU time = ' num2str(sols{i}.stats.cpu_time)],'Color',colors(col,:))
        legend('Location','northwest')
    end
end