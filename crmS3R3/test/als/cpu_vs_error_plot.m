function cpu_vs_error_plot(sols, varargin)
    % The first element of varargin should always be the variable of which we evaluate the global error!
    % If it is empty is assumed to be just q

    if nargin < 1
        error('Not enough arguments!')
    end
    
    size_var = max(size(varargin));
    if size_var > 0 && ~isempty(varargin{1})
        my_var = varargin{1};
    else
        my_var = 'q';
    end

    size_s = max(size(sols));
    colors = linspecer(size_s + 1);

    figure()
    hold on

    for i = 1:size_s
        if sols{i}.err.abs.(my_var) > 0
            text = [];
            for j = 1:size_var-1
                switch lower(varargin{j+1})
                case {'stab'}                    
                    if sols{i}.stab2 == 0
                        text = [text 'Original Solution'];
                    elseif sols{i}.stab_proj == 0
                        text = [text 'Baumgarte stabilization, $\gamma=$' num2str(sols{i}.a_baumgarte)];
                    else
                        text = [text 'Projection, $i_{\max}=$' num2str(sols{i}.imax)];
                    end
                case {'n'}
                    text = [text 'N=' num2str(sols{i}.n)];
                case {'step_control'}
                    if sols{i}.step_control == 1
                        text = [text 'variable time step size'];
                    else
                        text = [text 'fixed time step size'];
                    end
                otherwise
                    warning(['' varargin{j-1} ''' not yet implemented']);
                end
                if j < size_var-1
                    text = [text ', '];
                end                   
            end
            plot(sols{i}.err.rel.(my_var),sols{i}.stats.cpu_time,'o',...
                'MarkerSize', 5, 'MarkerEdgeColor', colors(i,:),...
                'MarkerFaceColor', colors(i,:), 'DisplayName',text);
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