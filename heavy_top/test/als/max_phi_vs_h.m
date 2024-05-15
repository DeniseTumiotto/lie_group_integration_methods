function max_phi_vs_h(sols, feature, side)

    % 'feature' : space, a_baumgarte
    % 'side' : 1 = TRUE, 0 = FALSE

    if nargin < 1
        error('Not enough arguments!')
    elseif nargin < 2
        % default value
        feature = 'a_baumgarte';
        side = 0;
    elseif nargin < 3
        side = 0;
    end

    if isempty(feature)
        feature = 'a_baumgarte';
    end

    size_s = max(size(sols));

    if ~isfield(sols{1},'h')
        sols = create_h(sols);
    end

    spaces = [];
    parameters = [];

    for i = 1:size_s
        if strcmp(feature,'a_baumgarte')
            if side
                if ~exist('n_fig','var')
                    figure()
                    n_fig = get(gcf,'Number');
                else
                    figure(n_fig)
                end
                k = find(spaces==sols{i}.problemset);
                if isempty(k)
                    spaces = [spaces sols{i}.problemset];
                    subplot(1,2,length(spaces))
                    plt{length(spaces)} = loglog(sols{i}.h, max(vecnorm(sols{i}.rslt.Phi)), '-', ...
                        'DisplayName', num2str(sols{i}.problemset), 'LineWidth', 1.5);
                else
                    subplot(1,2,k)
                    plt{k}.XData = [plt{k}.XData sols{i}.h];
                    plt{k}.YData = [plt{k}.YData max(vecnorm(sols{i}.rslt.Phi))];
                end
            else
                j = find(spaces==sols{i}.problemset);
                if isempty(j)
                    spaces = [spaces sols{i}.problemset];
                    figure(length(spaces))
                    hold on
                    k = find(parameters==sols{i}.a_baumgarte);
                    if isempty(k)
                        parameters = [parameters sols{i}.a_baumgarte];
                        plt{length(spaces),length(parameters)} = loglog(sols{i}.h, max(vecnorm(sols{i}.rslt.Phi)), 'o-', ...
                            'DisplayName', ['$\gamma$ = ' num2str(sols{i}.a_baumgarte)], 'LineWidth', 1.5);
                    else
                        plt{length(spaces),k} = loglog(sols{i}.h, max(vecnorm(sols{i}.rslt.Phi)), 'o-', ...
                            'DisplayName', ['$\gamma$ = ' num2str(sols{i}.a_baumgarte)], 'LineWidth', 1.5);
                    end
                else
                    figure(j)
                    k = find(parameters==sols{i}.a_baumgarte);
                    if isempty(k)
                        parameters = [parameters sols{i}.a_baumgarte];
                        plt{j,length(parameters)} = loglog(sols{i}.h, max(vecnorm(sols{i}.rslt.Phi)), 'o-', ...
                            'DisplayName', ['$\gamma$ = ' num2str(sols{i}.a_baumgarte)], 'LineWidth', 1.5);
                    else
                        if isempty(plt{j,k})
                            plt{j,k} = loglog(sols{i}.h, max(vecnorm(sols{i}.rslt.Phi)), 'o-', ...
                                'DisplayName', ['$\gamma$ = ' num2str(sols{i}.a_baumgarte)], 'LineWidth', 1.5);
                        else
                            plt{j,k}.XData = [plt{j,k}.XData sols{i}.h];
                            plt{j,k}.YData = [plt{j,k}.YData max(vecnorm(sols{i}.rslt.Phi))];
                        end
                    end
                end
            end
        elseif strcmp(feature,'space')
            disp('later')
        else
            disp('not implemented')
        end
        ax=gca;
        ax.FontSize = 14;
        ax.PlotBoxAspectRatio = [1 1 1];
        set(gca, 'YScale', 'log', 'XScale', 'log')
        switch sols{i}.liegroup
            case 1 % SO(3)
                disp('function not implemented for SO(3)')
            case 2 % S³
                disp('function not implemented for S^3')
            case 3 % SO(3) x R³
                space = '$SO(3)\times R^3$';
            case 4 % S³ x R³    
                space = '$S^3\times R^3$';
            case 5 % SE(3)
                space = '$SE(3)$';
            case 6 % S³ |x R³
                space = '$S^3\ \vert\times R^3$';
            case 7 % UDQ
        end
        title(['Max constraint residual in ', space],'FontSize',14,'Interpreter','latex')
        xlabel('Time step size (h)','FontSize',14,'Interpreter','latex')
        ylabel('$\max(\Vert\Phi(R,x)\Vert_2)$','FontSize',14,'Interpreter','latex')
        legend('Location','best','AutoUpdate','on','Interpreter','latex')
        grid on

    end
end