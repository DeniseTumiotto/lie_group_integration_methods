function cpu_vs_error_plot(sols, my_var)

    if nargin < 2
        my_var = 'q';
    end

    size_s = max(size(sols));
    colors = linspecer(min(4, size_s + 1)+1);

    figure()
    hold on

    for i = 1:size_s
        if sols{i}.err.abs.(my_var) > 0
            if sols{i}.stab2 == 0
                text = 'Original Solution';
            elseif sols{i}.stab_proj == 0
                text = ['Baumgarte stabilization, $\gamma=$' num2str(sols{i}.a_baumgarte)];
            else
                text = ['Projection, $i_{\max}=$' num2str(sols{i}.imax)];
            end
            plot(sols{i}.err.rel.(my_var),sols{i}.stats.cpu_time,'o',...
                'MarkerSize', 5, 'MarkerEdgeColor', colors(mod(i,size(colors,1))+1,:),...
                'MarkerFaceColor', colors(mod(i,size(colors,1))+1,:), 'DisplayName',text);
        end
    end

    ax=gca;
    ax.FontSize = 11; ax.PlotBoxAspectRatio = [1 1 1];
    set(gca, 'XScale', 'log')
    title('CPU time vs global error','FontSize',14,'Interpreter','latex')
    xlabel(['Global error in $' my_var '$'],'FontSize',14,'Interpreter','latex')
    ylabel('CPU time $(s)$','FontSize',14,'Interpreter','latex')
    legend('Location','best','AutoUpdate','on','Interpreter','latex','FontSize',14)
    grid on
    box on

end