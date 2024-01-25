function error_control_plot(sols, varargin)

if nargin < 1
    error('Not enough arguments!')
else
    for j = 2:2:nargin
        switch lower(varargin{j-1})
            case {'name'}
                the_name = varargin{j};
            case {'end time'}
                t_end = varargin{j};            
            case {'separate','alone'}
                separate = varargin{j};
            case {'combine','subplots'}
                combine = varargin{j};
            otherwise
                warning(['parameter ''' varargin{j-1} ''' not recognized']);
        end
    end
end

if ~exist('the_name','var')
    the_name = 'problem_name';
end
if ~exist('t_end','var')
    t_end = 1e7;
end
if ~exist('separate','var')
    separate = 0;
end
if ~exist('combine','var')
    combine = 1;
end

the_size = max(size(sols));
colors = linspecer(the_size+1);

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
        figure(i + the_size)
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