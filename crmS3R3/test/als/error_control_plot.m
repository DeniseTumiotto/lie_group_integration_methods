function error_control_plot(sols, combine, separate, t_end, the_name)

if nargin < 5
    the_name = 'problem_name';
    if nargin < 4
        t_end = 1e7;
        if nargin < 3
            separate = 0;
            if nargin < 2
                combine = 0;
            end
        end
    end
end

the_size = size(sols);
colors = linspecer(max(the_size)+1);

for i = 1:max(the_size)
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
    
    if ~combine
        figure(i)
    else
        figure(1)
    end
    if ~separate
        subplot(2,1,1)
    end
    plot(time_interval, sols{i}.rslt.local_err(2:j),...
        'LineWidth',1.5,'DisplayName',my_problem,'Color',colors(2*i,:))
    hold on
    title('\textbf{Variable} \textit{err} \textbf{for accepted time steps}',...
        'FontSize',14,'Interpreter','latex')

    xlabel('Time (t)',...
        'FontSize',14,'Interpreter','latex')

    ylabel('\textit{err}',...
        'FontSize',14,'Interpreter','latex')

    if separate
        legend()
    end

    if separate && ~combine
        figure(i + max(the_size))
    elseif separate && combine
        figure(2)
    else
        subplot(2,1,2)
    end
    plot(time_interval,sols{i}.h(1:j-1),...
        'LineWidth',1.5,'DisplayName',my_problem,'Color',colors(2*i,:))
    hold on
    title('\textbf{Time step size}',...
        'FontSize',14,'Interpreter','latex')

    xlabel('Time (t)',...
        'FontSize',14,'Interpreter','latex')

    ylabel('$h_{\textrm{new}}$',...
        'FontSize',14,'Interpreter','latex')

    legend()

end

end