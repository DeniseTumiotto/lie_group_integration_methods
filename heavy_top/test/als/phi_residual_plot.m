function phi_residual_plot(sols, varargin)

    % varargin could take as values:
    % 'side'
    % 'overlap'
    % 'feature' : space, a_baumgarte, h
    % after the keyword write the value: 1 for TRUE, 0 for FALSE

    % default values
    side    = 0;
    overlap = 0;
    feature = 'a_baumgarte';

    if nargin < 1
        error('Not enough arguments!')
    elseif nargin > 1
        for j = 2:2:nargin
            switch lower(varargin{j-1})
                case {'side'}
                    side = varargin{j};
                case{'overlap'}
                    overlap = varargin{j};
                case{'feature'}
                    feature = varargin{j};
                otherwise
                    warning(['parameter ''' varargin{j-1} ''' not recognized']);
            end
        end
    end

    size_s     = max(size(sols));
    parameters = [];
    spaces     = [];
    index_cntr = {};

    % if one wants the plot for each time step size but it is not evaluated
    if strcmp(feature,'h') && ~isfield(sols{1},'h')
        sols = create_h(sols);
    elseif strcmp(feature,'space')
        feature = 'problemset';
    end

    for i = 1:size_s
        l = find(spaces==sols{i}.problemset);
        if isempty(l)
            spaces = [spaces sols{i}.problemset];
            index_cntr{end+1} = i;
        else
            index_cntr{l} = [index_cntr{l} i];
        end
    end

    l_spaces = length(spaces);

    for i = 1:l_spaces
        if side
            if overlap
                if l_spaces > 2
                    figure(1)
                    if sols{index_cntr{i}(1)}.liegroup == 3 || sols{index_cntr{i}(1)}.liegroup == 5
                        subplot(1,2,1)
                        hold on
                    else
                        subplot(1,2,2)
                    end
                else
                    disp('not possible!')
                    break
                end
            else
                figure(1)
                subplot(l_spaces/2,2,i)
                hold on
            end
        else
            if overlap
                if sols{index_cntr{i}(1)}.liegroup == 3 || sols{index_cntr{i}(1)}.liegroup == 5
                    figure(1)
                    hold on
                else
                    figure(2)
                    hold on
                end
            else
                figure()
                hold on
            end
        end
        for j = index_cntr{i}
            plot(sols{j}.rslt.t, vecnorm(sols{j}.rslt.Phi), ...
                'DisplayName', [num2str(spaces(i)), feature '=' num2str(sols{i}.(feature))], 'LineWidth', 1.5)
        end
    end
    legend('Location','best','AutoUpdate','on','Interpreter','latex')
end