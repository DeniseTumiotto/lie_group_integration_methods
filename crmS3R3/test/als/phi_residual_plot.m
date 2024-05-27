function phi_residual_plot(sols)

    size_s = max(size(sols));
    colors = linspecer(min(4, size_s + 1)+1);
    style  = ['-',':','--','-.'];

    figure()
    hold on
    for i = 1:size_s
        plot(sols{i}.rslt.t, vecnorm(sols{i}.rslt.Phi), ...
            'Color', colors(mod(i,size(colors,1))+1,:), 'LineWidth', 1.5, ...
            'LineStyle', style(mod(i,4)+1))
    end
    ax=gca;
    ax.FontSize = 14;
    ax.PlotBoxAspectRatio = [1 1 1];
    ax.YScale = 'log';
    title(['Constraint residual vs time'],'FontSize',14,'Interpreter','latex')
    xlabel('Time $(t)$','FontSize',14,'Interpreter','latex')
    ylabel('$\Vert\Phi(R,x)\Vert_2$','FontSize',14,'Interpreter','latex')
    grid on
    box on
end