function phi_residual_plot(sols,name)

    % name: C
    % name: h
    % name: a_baumgarte

    if nargin < 2
        name = 'a_baumgarte';
    end

    size_s = max(size(sols));
    colors = linspecer(min(4, size_s + 1)+1);
    style  = {'-', ':', '--', '-.'};

    if ~isfield(sols{1},'h')
        sols = create_h(sols);
    end
    
    figure()
    hold on
    for i = 1:size_s
        if strcmp(name, 'a_baumgarte')
            text = ['$\gamma =$' num2str(sols{i}.(name))];
        elseif strcmp(name, 'C')
            text = ['$C=$' num2str(sols{i}.a_baumgarte*sols{i}.h(1))];
        elseif strcmp(name, 'h')
            text = ['$h=$' num2str(sols{i}.h)];
        else
            text = '';
        end

        plot(sols{i}.rslt.t, vecnorm(sols{i}.rslt.Phi), ...
            'DisplayName', text, 'Color', colors(mod(i,size(colors,1))+1,:), ...
            'LineWidth', 1.5, 'LineStyle', style{mod(i,4)+1})
    end
    ax=gca;
    ax.FontSize = 14;
    ax.PlotBoxAspectRatio = [1 1 1];
    ax.YScale = 'log';
    title('Constraint residual vs time','FontSize',14,'Interpreter','latex')
    xlabel('Time $(t)$','FontSize',14,'Interpreter','latex')
    ylabel('$\Vert\Phi(R,x)\Vert_2$','FontSize',14,'Interpreter','latex')
    legend('Interpreter','latex')
    grid on
    box on
end