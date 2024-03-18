function phi_residual_plot(sols, side)

size_s = max(size(sols));

if nargin < 2
    side = 0;
end

for i = 1:size_s
    switch sols{i}.liegroup
        case 1 % SO(3)
        case 2 % S³
        case 3 % SO(3) x R³
            if side
                if ~exist('n_fig','var')
                    figure()
                    n_fig = get(gcf,'Number');
                    subplot(1,2,1)
                else
                    figure(n_fig)
                    subplot(1,2,1)
                end
            else
                if ~exist('n_so3r3','var')
                    figure()
                    n_so3r3 = get(gcf,'Number');
                else
                    figure(n_so3r3)
                end
            end
            space = '$SO(3)\times R^3$';
        case 4 % S³ x R³
            if ~exist('n_s3r3','var')
                figure()
                n_s3r3 = get(gcf,'Number');
            else
                figure(n_s3r3)
            end
            space = '$S^3\times R^3$';
        case 5 % SE(3)
            if side
                if ~exist('n_fig','var')
                    figure()
                    n_fig = get(gcf,'Number');
                    subplot(1,2,2)
                else
                    figure(n_fig)
                    subplot(1,2,2)
                end
            else
                if ~exist('n_se3','var')
                    figure()
                    n_se3 = get(gcf,'Number');
                else
                    figure(n_se3)
                end
            end
            space = '$SE(3)$';
        case 6 % S³ |x R³            
            if ~exist('n_s3semi','var')
                figure()
                n_s3semi = get(gcf,'Number');
            else
                figure(n_s3semi)
            end
            space = '$SO(3)\ \vert\times R^3$';
        case 7 % UDQ
    end
    hold on
    ax=gca;
    ax.FontSize = 14;
    ax.PlotBoxAspectRatio = [1 1 1];
    title(['Constraint residual in ', space],'FontSize',14,'Interpreter','latex')
    xlabel('Time (t)','FontSize',14,'Interpreter','latex')
    ylabel('$\Vert\Phi(R,x)\Vert_2$','FontSize',14,'Interpreter','latex')
    legend('Location','best','AutoUpdate','on')
    grid on
    plot(sols{i}.rslt.t,vecnorm(sols{i}.rslt.Phi),'DisplayName',['a = ' num2str(sols{i}.a_baumgarte)],'LineWidth',1.5)
end

end