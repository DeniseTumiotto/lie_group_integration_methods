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

    size_s = max(size(sols));

    % if one wants the plot for each time step size but it is not evaluated
    if strcmp(feature,'h') && ~isfield(sols{1},'h')
        sols = create_h(sols);
    elseif strcmp(feature,'space')
        feature = 'problemset';
    end

    for i = 1:size_s
        switch sols{i}.liegroup
            case 1 % SO(3)
                disp('function not implemented for SO(3)')
            case 2 % S³
                disp('function not implemented for S^3')
            case 3 % SO(3) x R³
                if side
                    if ~exist('n_fig','var')
                        figure()
                        n_fig = get(gcf,'Number');
                        if overlap
                            subplot(1,2,1)
                            line_name = '$SO(3)\times R^3$, ';
                        else
                            subplot(2,2,1)
                            line_name = '';
                        end
                    else
                        figure(n_fig)
                        if overlap
                            subplot(1,2,1)
                            line_name = '$SO(3)\times R^3$, ';
                        else
                            subplot(2,2,1)
                            line_name = '';
                        end
                    end
                else
                    if ~exist('n_so3r3','var') && ~exist('n_so3','var')
                        figure()
                        if overlap
                            line_name = '$SO(3)\times R^3$, ';
                            n_so3 = get(gcf,'Number');
                        else
                            line_name = '';
                            n_so3r3 = get(gcf,'Number');
                        end
                    else
                        if overlap
                            figure(n_so3)
                            line_name = '$SO(3)\times R^3$, ';
                        else
                            figure(n_so3r3)
                            line_name = '';
                        end
                    end
                end
                space = '$SO(3)\times R^3$';
            case 4 % S³ x R³
                if side
                    if ~exist('n_fig','var')
                        figure()
                        n_fig = get(gcf,'Number');
                        if overlap
                            subplot(1,2,2)
                            line_name = '$S^3\times R^3$, ';
                        else
                            subplot(2,2,3)
                            line_name = '';
                        end
                    else
                        figure(n_fig)
                        if overlap
                            subplot(1,2,2)
                            line_name = '$S^3\times R^3$, ';
                        else
                            subplot(2,2,3)
                            line_name = '';
                        end
                    end
                else
                    if ~exist('n_s3r3','var') && ~exist('n_s3','var')
                        figure()
                        if overlap
                            line_name = '$S^3\times R^3$, ';
                            n_s3 = get(gcf,'Number');
                        else
                            line_name = '';
                            n_s3r3 = get(gcf,'Number');
                        end
                    else
                        if overlap
                            figure(n_s3)
                            line_name = '$S^3\times R^3$, ';
                        else
                            figure(n_s3r3)
                            line_name = '';
                        end
                    end
                end            
                space = '$S^3\times R^3$';
            case 5 % SE(3)
                if side
                    if ~exist('n_fig','var')
                        figure()
                        n_fig = get(gcf,'Number');
                        if overlap
                            subplot(1,2,1)
                            line_name = '$SE3$, ';
                        else
                            subplot(2,2,2)
                            line_name = '';
                        end
                    else
                        figure(n_fig)
                        if overlap
                            subplot(1,2,1)
                            line_name = '$SE3$, ';
                        else
                            subplot(2,2,2)
                            line_name = '';
                        end
                    end
                else
                    if ~exist('n_se3','var') && ~exist('n_so3','var')
                        figure()
                        if overlap
                            line_name = '$SE3$, ';
                            n_so3 = get(gcf,'Number');
                        else
                            line_name = '';
                            n_se3 = get(gcf,'Number');
                        end
                    else
                        if overlap
                            figure(n_so3)
                            line_name = '$SE3$, ';
                        else
                            figure(n_se3)
                            line_name = '';
                        end
                    end
                end
                space = '$SE(3)$';
            case 6 % S³ |x R³
                if side
                    if ~exist('n_fig','var')
                        figure()
                        n_fig = get(gcf,'Number');
                        if overlap
                            subplot(1,2,2)
                            line_name = '$S^3\vert\times R^3$, ';
                        else
                            subplot(2,2,4)
                            line_name = '';
                        end
                    else
                        figure(n_fig)
                        if overlap
                            subplot(1,2,2)
                            line_name = '$S^3\vert\times R^3$, ';
                        else
                            subplot(2,2,4)
                            line_name = '';
                        end
                    end
                else
                    if ~exist('n_s3sdr3','var') && ~exist('n_s3','var')
                        figure()
                        if overlap
                            line_name = '$S^3\vert\times R^3$, ';
                            n_s3 = get(gcf,'Number');
                        else
                            line_name = '';
                            n_s3sdr3 = get(gcf,'Number');
                        end
                    else
                        if overlap
                            figure(n_s3)
                            line_name = '$S^3\vert\times R^3$, ';
                        else
                            figure(n_s3sdr3)
                            line_name = '';
                        end
                    end
                end
                space = '$S^3\ \vert\times R^3$';
            case 7 % UDQ
        end
        hold on
        ax=gca;
        ax.FontSize = 14;
        ax.PlotBoxAspectRatio = [1 1 1];
        if overlap
            title(['Constraint residual vs time (t)'],'FontSize',14,'Interpreter','latex')
        else
            title(['Constraint residual in ', space],'FontSize',14,'Interpreter','latex')
        end
        xlabel('Time (t)','FontSize',14,'Interpreter','latex')
        ylabel('$\Vert\Phi(R,x)\Vert_2$','FontSize',14,'Interpreter','latex')
        grid on
        plot(sols{i}.rslt.t, vecnorm(sols{i}.rslt.Phi), 'DisplayName', [line_name feature '=' num2str(sols{i}.(feature))], 'LineWidth', 1.5)
        legend('Location','best','AutoUpdate','on')
    end

end