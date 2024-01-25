function orders_plot(sols, varargin)

if nargin < 1
    error('Not enough arguments!')
else
    for j = 2:2:nargin
        switch lower(varargin{j-1})
            case {'variable'}
                my_var = varargin{j};
            otherwise
                warning(['parameter ''' varargin{j-1} ''' not recognized']);
        end
    end
end
if ~exist('my_var','var')
    my_var = 'all';
end

the_size = max(size(sols));
colors = linspecer(2 * 4);
the_style = {'-';'-.';':'};

for i = 1:the_size
    if exist('n_figure','var')
        figure(n_figure)
    else
        figure()
        hold on
        n_figure = get(gcf,'Number');
        ax=gca;
        ax.FontSize = 14;
        set(gca, 'XScale', 'log', 'YScale', 'log')
        if length(my_var) == 1
            title(['Order of convergence, ' my_var],'FontSize',14)
        else
            title('Order of convergence','FontSize',14)
        end
        xlabel('Time step size (h)','FontSize',14)
        ylabel('Global error','FontSize',14)
        legend('Location','best','AutoUpdate','on')
        grid on
    end
    switch my_var
        case 'all'
                switch sols{i}.problem_name
                    case 'roll-up'
                        if exist('pl_rll','var') && ~isempty(pl_rll{sols{i}.order,1})
                            pl_rll{sols{i}.order,1}.XData = [pl_rll{sols{i}.order,1}.XData,sols{i}.h];
                            pl_rll{sols{i}.order,1}.YData = [pl_rll{sols{i}.order,1}.YData,sols{i}.err.abs.q];
                            pl_rll{sols{i}.order,2}.XData = [pl_rll{sols{i}.order,2}.XData,sols{i}.h];
                            pl_rll{sols{i}.order,2}.YData = [pl_rll{sols{i}.order,2}.YData,sols{i}.err.abs.v];
                            pl_rll{sols{i}.order,3}.XData = [pl_rll{sols{i}.order,3}.XData,sols{i}.h];
                            pl_rll{sols{i}.order,3}.YData = [pl_rll{sols{i}.order,3}.YData,sols{i}.err.abs.l];
                        else
                            pl_rll = {[], [], [];
                                      [], [], [];
                                      [], [], [];
                                      [], [], [];
                                      [], [], []};
                            pl_rll{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.q,...
                                'Color',colors(- 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle',the_style{1}, ...
                                'DisplayName',['q, order: ' num2str(sols{i}.order) ', Roll-up']);
                            pl_rll{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.v,...
                                'Color',colors(- 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{2}, ...
                                'DisplayName',['v, order: ' num2str(sols{i}.order) ', Roll-up']);
                            pl_rll{sols{i}.order,3} = loglog(sols{i}.h, sols{i}.err.abs.l,...
                                'Color',colors(- 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{3}, ...
                                'DisplayName',['l, order: ' num2str(sols{i}.order) ', Roll-up']);
                        end
                    case 'flying_spaghetti'
                        if exist('pl_fly','var') && ~isempty(pl_fly{sols{i}.order,1})
                            pl_fly{sols{i}.order,1}.XData = [pl_fly{sols{i}.order,1}.XData,sols{i}.h];
                            pl_fly{sols{i}.order,1}.YData = [pl_fly{sols{i}.order,1}.YData,sols{i}.err.abs.q];
                            pl_fly{sols{i}.order,2}.XData = [pl_fly{sols{i}.order,2}.XData,sols{i}.h];
                            pl_fly{sols{i}.order,2}.YData = [pl_fly{sols{i}.order,2}.YData,sols{i}.err.abs.v];
                            pl_fly{sols{i}.order,3}.XData = [pl_fly{sols{i}.order,3}.XData,sols{i}.h];
                            pl_fly{sols{i}.order,3}.YData = [pl_fly{sols{i}.order,3}.YData,sols{i}.err.abs.l];
                        else
                            pl_fly = {[], [], [];
                                      [], [], [];
                                      [], [], [];
                                      [], [], [];
                                      [], [], []};
                            pl_fly{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.q,...
                                'Color',colors(4 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle',the_style{1}, ...
                                'DisplayName',['q, order: ' num2str(sols{i}.order) ', Flying spaghetti']);
                            pl_fly{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.v,...
                                'Color',colors(4 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{2}, ...
                                'DisplayName',['v, order: ' num2str(sols{i}.order) ', Flying spaghetti']);
                            pl_fly{sols{i}.order,3} = loglog(sols{i}.h, sols{i}.err.abs.l,...
                                'Color',colors(4 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{3}, ...
                                'DisplayName',['l, order: ' num2str(sols{i}.order) ', Flying spaghetti']);
                        end
                    otherwise
                        if exist('pl_oth','var') && ~isempty(pl_oth{sols{i}.order,1})
                            pl_oth{sols{i}.order,1}.XData = [pl_oth{sols{i}.order,1}.XData,sols{i}.h];
                            pl_oth{sols{i}.order,1}.YData = [pl_oth{sols{i}.order,1}.YData,sols{i}.err.abs.q];
                            pl_oth{sols{i}.order,2}.XData = [pl_oth{sols{i}.order,2}.XData,sols{i}.h];
                            pl_oth{sols{i}.order,2}.YData = [pl_oth{sols{i}.order,2}.YData,sols{i}.err.abs.v];
                            pl_oth{sols{i}.order,3}.XData = [pl_oth{sols{i}.order,3}.XData,sols{i}.h];
                            pl_oth{sols{i}.order,3}.YData = [pl_oth{sols{i}.order,3}.YData,sols{i}.err.abs.l];
                        else
                            pl_oth = {[], [], [];
                                      [], [], [];
                                      [], [], [];
                                      [], [], [];
                                      [], [], []};
                            pl_oth{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.q,...
                                'Color',colors(8 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle',the_style{1}, ...
                                'DisplayName',['q, order: ' num2str(sols{i}.order) ', Other']);
                            pl_oth{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.v,...
                                'Color',colors(8 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{2}, ...
                                'DisplayName',['v, order: ' num2str(sols{i}.order) ', Other']);
                            pl_oth{sols{i}.order,3} = loglog(sols{i}.h, sols{i}.err.abs.l,...
                                'Color',colors(8 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{3}, ...
                                'DisplayName',['l, order: ' num2str(sols{i}.order) ', Other']);
                        end
                end
        case 'q v'
            switch sols{i}.problem_name
                case 'roll-up'
                    if exist('pl_rll','var') && ~isempty(pl_rll{sols{i}.order,1})
                        pl_rll{sols{i}.order,1}.XData = [pl_rll{sols{i}.order,1}.XData,sols{i}.h];
                        pl_rll{sols{i}.order,1}.YData = [pl_rll{sols{i}.order,1}.YData,sols{i}.err.abs.q];
                        pl_rll{sols{i}.order,2}.XData = [pl_rll{sols{i}.order,2}.XData,sols{i}.h];
                        pl_rll{sols{i}.order,2}.YData = [pl_rll{sols{i}.order,2}.YData,sols{i}.err.abs.v];
                    else
                        pl_rll = {[], [];
                                    [], [];
                                    [], [];
                                    [], [];
                                    [], []};
                        pl_rll{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.q,...
                            'Color',colors(- 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle',the_style{1}, ...
                            'DisplayName',['q, order: ' num2str(sols{i}.order) ', Roll-up']);
                        pl_rll{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.v,...
                            'Color',colors(- 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{2}, ...
                            'DisplayName',['v, order: ' num2str(sols{i}.order) ', Roll-up']);
                    end
                case 'flying_spaghetti'
                    if exist('pl_fly','var') && ~isempty(pl_fly{sols{i}.order,1})
                        pl_fly{sols{i}.order,1}.XData = [pl_fly{sols{i}.order,1}.XData,sols{i}.h];
                        pl_fly{sols{i}.order,1}.YData = [pl_fly{sols{i}.order,1}.YData,sols{i}.err.abs.q];
                        pl_fly{sols{i}.order,2}.XData = [pl_fly{sols{i}.order,2}.XData,sols{i}.h];
                        pl_fly{sols{i}.order,2}.YData = [pl_fly{sols{i}.order,2}.YData,sols{i}.err.abs.v];
                    else
                        pl_fly = {[], [];
                                    [], [];
                                    [], [];
                                    [], [];
                                    [], []};
                        pl_fly{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.q,...
                            'Color',colors(4 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle',the_style{1}, ...
                            'DisplayName',['q, order: ' num2str(sols{i}.order) ', Flying spaghetti']);
                        pl_fly{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.v,...
                            'Color',colors(4 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{2}, ...
                            'DisplayName',['v, order: ' num2str(sols{i}.order) ', Flying spaghetti']);
                    end
                otherwise
                    if exist('pl_oth','var') && ~isempty(pl_oth{sols{i}.order,1})
                        pl_oth{sols{i}.order,1}.XData = [pl_oth{sols{i}.order,1}.XData,sols{i}.h];
                        pl_oth{sols{i}.order,1}.YData = [pl_oth{sols{i}.order,1}.YData,sols{i}.err.abs.q];
                        pl_oth{sols{i}.order,2}.XData = [pl_oth{sols{i}.order,2}.XData,sols{i}.h];
                        pl_oth{sols{i}.order,2}.YData = [pl_oth{sols{i}.order,2}.YData,sols{i}.err.abs.v];
                    else
                        pl_oth = {[], [];
                                    [], [];
                                    [], [];
                                    [], [];
                                    [], []};
                        pl_oth{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.q,...
                            'Color',colors(8 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle',the_style{1}, ...
                            'DisplayName',['q, order: ' num2str(sols{i}.order) ', Other']);
                        pl_oth{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.v,...
                            'Color',colors(8 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{2}, ...
                            'DisplayName',['v, order: ' num2str(sols{i}.order) ', Other']);
                    end
            end
        case 'v l'
            switch sols{i}.problem_name
                case 'roll-up'
                    if exist('pl_rll','var') && ~isempty(pl_rll{sols{i}.order,1})
                        pl_rll{sols{i}.order,1}.XData = [pl_rll{sols{i}.order,1}.XData,sols{i}.h];
                        pl_rll{sols{i}.order,1}.YData = [pl_rll{sols{i}.order,1}.YData,sols{i}.err.abs.v];
                        pl_rll{sols{i}.order,2}.XData = [pl_rll{sols{i}.order,2}.XData,sols{i}.h];
                        pl_rll{sols{i}.order,2}.YData = [pl_rll{sols{i}.order,2}.YData,sols{i}.err.abs.l];
                    else
                        pl_rll = {[], [];
                                    [], [];
                                    [], [];
                                    [], [];
                                    [], []};
                        pl_rll{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.v,...
                            'Color',colors(- 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{2}, ...
                            'DisplayName',['v, order: ' num2str(sols{i}.order) ', Roll-up']);
                        pl_rll{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.l,...
                            'Color',colors(- 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{3}, ...
                            'DisplayName',['l, order: ' num2str(sols{i}.order) ', Roll-up']);
                    end
                case 'flying_spaghetti'
                    if exist('pl_fly','var') && ~isempty(pl_fly{sols{i}.order,1})
                        pl_fly{sols{i}.order,1}.XData = [pl_fly{sols{i}.order,1}.XData,sols{i}.h];
                        pl_fly{sols{i}.order,1}.YData = [pl_fly{sols{i}.order,1}.YData,sols{i}.err.abs.v];
                        pl_fly{sols{i}.order,2}.XData = [pl_fly{sols{i}.order,2}.XData,sols{i}.h];
                        pl_fly{sols{i}.order,2}.YData = [pl_fly{sols{i}.order,2}.YData,sols{i}.err.abs.l];
                    else
                        pl_fly = {[], [];
                                    [], [];
                                    [], [];
                                    [], [];
                                    [], []};
                        pl_fly{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.v,...
                            'Color',colors(4 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{2}, ...
                            'DisplayName',['v, order: ' num2str(sols{i}.order) ', Flying spaghetti']);
                        pl_fly{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.l,...
                            'Color',colors(4 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{3}, ...
                            'DisplayName',['l, order: ' num2str(sols{i}.order) ', Flying spaghetti']);
                    end
                otherwise
                    if exist('pl_oth','var') && ~isempty(pl_oth{sols{i}.order,1})
                        pl_oth{sols{i}.order,1}.XData = [pl_oth{sols{i}.order,1}.XData,sols{i}.h];
                        pl_oth{sols{i}.order,1}.YData = [pl_oth{sols{i}.order,1}.YData,sols{i}.err.abs.v];
                        pl_oth{sols{i}.order,2}.XData = [pl_oth{sols{i}.order,2}.XData,sols{i}.h];
                        pl_oth{sols{i}.order,2}.YData = [pl_oth{sols{i}.order,2}.YData,sols{i}.err.abs.l];
                    else
                        pl_oth = {[], [];
                                    [], [];
                                    [], [];
                                    [], [];
                                    [], []};
                        pl_oth{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.v,...
                            'Color',colors(8 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{2}, ...
                            'DisplayName',['v, order: ' num2str(sols{i}.order) ', Other']);
                        pl_oth{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.l,...
                            'Color',colors(8 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{3}, ...
                            'DisplayName',['l, order: ' num2str(sols{i}.order) ', Other']);
                    end
            end
        case 'q l'
            switch sols{i}.problem_name
                case 'roll-up'
                    if exist('pl_rll','var') && ~isempty(pl_rll{sols{i}.order,1})
                        pl_rll{sols{i}.order,1}.XData = [pl_rll{sols{i}.order,1}.XData,sols{i}.h];
                        pl_rll{sols{i}.order,1}.YData = [pl_rll{sols{i}.order,1}.YData,sols{i}.err.abs.q];
                        pl_rll{sols{i}.order,2}.XData = [pl_rll{sols{i}.order,2}.XData,sols{i}.h];
                        pl_rll{sols{i}.order,2}.YData = [pl_rll{sols{i}.order,2}.YData,sols{i}.err.abs.l];
                    else
                        pl_rll = {[], [];
                                    [], [];
                                    [], [];
                                    [], [];
                                    [], []};
                        pl_rll{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.q,...
                            'Color',colors(- 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle',the_style{1}, ...
                            'DisplayName',['q, order: ' num2str(sols{i}.order) ', Roll-up']);
                        pl_rll{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.l,...
                            'Color',colors(- 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{3}, ...
                            'DisplayName',['l, order: ' num2str(sols{i}.order) ', Roll-up']);
                    end
                case 'flying_spaghetti'
                    if exist('pl_fly','var') && ~isempty(pl_fly{sols{i}.order,1})
                        pl_fly{sols{i}.order,1}.XData = [pl_fly{sols{i}.order,1}.XData,sols{i}.h];
                        pl_fly{sols{i}.order,1}.YData = [pl_fly{sols{i}.order,1}.YData,sols{i}.err.abs.q];
                        pl_fly{sols{i}.order,2}.XData = [pl_fly{sols{i}.order,2}.XData,sols{i}.h];
                        pl_fly{sols{i}.order,2}.YData = [pl_fly{sols{i}.order,2}.YData,sols{i}.err.abs.l];
                    else
                        pl_fly = {[], [];
                                    [], [];
                                    [], [];
                                    [], [];
                                    [], []};
                        pl_fly{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.q,...
                            'Color',colors(4 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle',the_style{1}, ...
                            'DisplayName',['q, order: ' num2str(sols{i}.order) ', Flying spaghetti']);
                        pl_fly{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.l,...
                            'Color',colors(4 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{3}, ...
                            'DisplayName',['l, order: ' num2str(sols{i}.order) ', Flying spaghetti']);
                    end
                otherwise
                    if exist('pl_oth','var') && ~isempty(pl_oth{sols{i}.order,1})
                        pl_oth{sols{i}.order,1}.XData = [pl_oth{sols{i}.order,1}.XData,sols{i}.h];
                        pl_oth{sols{i}.order,1}.YData = [pl_oth{sols{i}.order,1}.YData,sols{i}.err.abs.q];
                        pl_oth{sols{i}.order,2}.XData = [pl_oth{sols{i}.order,2}.XData,sols{i}.h];
                        pl_oth{sols{i}.order,2}.YData = [pl_oth{sols{i}.order,2}.YData,sols{i}.err.abs.l];
                    else
                        pl_oth = {[], [];
                                    [], [];
                                    [], [];
                                    [], [];
                                    [], []};
                        pl_oth{sols{i}.order,1} = loglog(sols{i}.h, sols{i}.err.abs.q,...
                            'Color',colors(8 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle',the_style{1}, ...
                            'DisplayName',['q, order: ' num2str(sols{i}.order) ', Other']);
                        pl_oth{sols{i}.order,2} = loglog(sols{i}.h, sols{i}.err.abs.l,...
                            'Color',colors(8 - 1 + sols{i}.order,:),'LineWidth',1.5,'LineStyle', the_style{3}, ...
                            'DisplayName',['l, order: ' num2str(sols{i}.order) ', Other']);
                    end
            end
        otherwise
            if length(my_var) == 1
                switch sols{i}.problem_name
                    case 'roll-up'
                        if exist('pl_rll','var') && ~isempty(pl_rll{sols{i}.order})
                            pl_rll{sols{i}.order}.XData = [pl_rll{sols{i}.order}.XData,sols{i}.h];
                            pl_rll{sols{i}.order}.YData = [pl_rll{sols{i}.order}.YData,sols{i}.err.abs.(my_var)];
                        else
                            pl_rll = {[], [], [], [], []};
                            pl_rll{sols{i}.order} = loglog(sols{i}.h,sols{i}.err.abs.(my_var),...
                                'Color',colors(- 1 + sols{i}.order,:),'LineWidth',1.5,'DisplayName',['order: ' num2str(sols{i}.order) ', Roll-up' ]);
                        end
                    case 'flying_spaghetti'
                        if exist('pl_fly','var') && ~isempty(pl_fly{sols{i}.order})
                            pl_fly{sols{i}.order}.XData = [pl_fly{sols{i}.order}.XData,sols{i}.h];
                            pl_fly{sols{i}.order}.YData = [pl_fly{sols{i}.order}.YData,sols{i}.err.abs.(my_var)];
                        else
                            pl_fly = {[], [], [], [], []};
                            pl_fly{sols{i}.order} = loglog(sols{i}.h,sols{i}.err.abs.(my_var),...
                                'Color',colors(4 - 1 + sols{i}.order,:),'LineWidth',1.5,'DisplayName',['order: ' num2str(sols{i}.order) ', Flying spaghetti']);
                        end
                    otherwise
                        if exist('pl_oth','var') && ~isempty(pl_oth{sols{i}.order})
                            pl_oth{sols{i}.order}.XData = [pl_oth{sols{i}.order}.XData,sols{i}.h];
                            pl_oth{sols{i}.order}.YData = [pl_oth{sols{i}.order}.YData,sols{i}.err.abs.(my_var)];
                        else
                            pl_oth = {[], [], [], [], []};
                            pl_oth{sols{i}.order} = loglog(sols{i}.h,sols{i}.err.abs.(my_var),...
                                'Color',colors(8 - 1 + sols{i}.order,:),'LineWidth',1.5,'DisplayName',['order: ' num2str(sols{i}.order) ', Other tests' ]);
                        end
                end
            else
                warning(['variable ''' my_var ''' not recognized!'])
            end
    end
end

end