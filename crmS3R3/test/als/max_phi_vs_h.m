function max_phi_vs_h(sols,name)

    % name: C
    % name: a_baumgarte

    if nargin < 2
        name = 'a_baumgarte';
    end

    size_s = max(size(sols));

    if ~isfield(sols{1},'h')
        sols = create_h(sols);
    end

    colors = linspecer(5);
    parameter = [];
    tol = 1e-4;

    figure()
    hold on

    for i = 1:size_s

        if strcmp(name, 'a_baumgarte')
            the_par = sols{i}.a_baumgarte;
            k = find(parameter==the_par);
            text = ['$\gamma =$' num2str(sols{i}.(name))];
        elseif strcmp(name, 'C')
            the_par = sols{i}.a_baumgarte*sols{i}.h;
            k = find(parameter<the_par+tol & parameter>the_par-tol);
            text = ['$C=$' num2str(sols{i}.a_baumgarte*sols{i}.h)];
        else
            text = '';
            k = [];
        end

        if isempty(k)
            parameter = [parameter the_par];
            plt{length(parameter)} = loglog(sols{i}.h, max(vecnorm(sols{i}.rslt.Phi)), 'o-', ...
                'DisplayName', text, 'LineWidth', 1.5, 'Color', colors(mod(i,k)+1,:));
        else
            plt{k}.XData = [plt{k}.XData sols{i}.h];
            plt{k}.YData = [plt{k}.YData max(vecnorm(sols{i}.rslt.Phi))];
        end
    end

    ax=gca;
    ax.FontSize = 14;
    ax.PlotBoxAspectRatio = [1 1 1];
    set(gca, 'YScale', 'log', 'XScale', 'log')

    title(['{\bf{Max constraint residual vs time step size}} $h$'],'FontSize',14,'Interpreter','latex')
    xlabel('Time step size ($h$)','FontSize',14,'Interpreter','latex')
    ylabel('$\max(\Vert{\bf{\Phi}}(q)\Vert_2)$','FontSize',14,'Interpreter','latex')
    legend('Location','best','AutoUpdate','on','Interpreter','latex')
    grid on
end