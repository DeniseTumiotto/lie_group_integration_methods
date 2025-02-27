function plot_confront(sols,interp)
figure()
hold on;

n = max(size(sols));
color = linspecer(n+1);

for it=1:n
    if interp
        [~, xs] = q_to_ps_xs(sols{it}.rslt.q_interp(:), sols{end},1);
    else
        [~, xs] = q_to_ps_xs(sols{it}.rslt.q(:,end), sols{it},0);
    end

   % Plot
   plot3(xs(1,:), xs(2,:), xs(3,:), '-o', 'MarkerSize',5, 'LineWidth',3, ...
         'DisplayName', num2str(sols{it}.n),'MarkerFaceColor', color(end-(it-1),:),...
         'Color', color(end-(it-1),:));
end
ax=gca;
ax.FontSize = 16;
ax.DataAspectRatio = [1,1,1];
box on
if strcmp(sols{1}.problem_name,'roll-up')
    view([0 -1 0])
    ax.XLim = [-4, 10];
    ax.ZLim = [-6, 2];
elseif strcmp(sols{1}.problem_name,'flying_spaghetti')
    view([22.103969891487075,21.233658241582319])
else
    view(2)
end

legend('FontSize', 16, 'Location','southeast');
xlabel('$x$', 'FontSize', 24, 'Interpreter', 'latex');
ylabel('$y$', 'FontSize', 24, 'Interpreter', 'latex');
zlabel('$z$', 'FontSize', 24, 'Interpreter', 'latex');

grid off;
grid on;

end

function [ps, xs] = q_to_ps_xs(q,sol,interp)

ps = zeros(4,sol.n+1);
xs = zeros(3,sol.n+1);

    if interp
        tmp = reshape(q,[7,sol.n+1]);
        ps = tmp(1:4,:);
        xs = tmp(5:7,:);
    else
if sol.output_s_at == 0

   if sol.fixed_x0 == 0
      xs(:,1) = q(5:7);
      q_start = 8;
   else
      xs(:,1) = sol.fixed_x0_position;
      q_start = 5;
   end
   if sol.fixed_p0 == 1
      ps(:,1) = sol.fixed_p0_orientation;
      q_start = q_start - 4;
   else
      ps(:,1) = q(1:4);
   end
   if sol.fixed_xn == 0
      xs(:,end) = q(end-2:end);
      ps(:,end) = q(end-6:end-3);
      q_end = length(q) - 7;
   else
      xs(:,end) = sol.fixed_xn_position;
      ps(:,end) = q(end-3:end);
      q_end = length(q) - 4;
   end
   if sol.fixed_pn == 1
      ps(:,end) = sol.fixed_pn_orientation;
      q_end = q_end - 4;
   end

   tmp = reshape(q(q_start:q_start-1 + 7*(sol.n-1)),[7,sol.n-1]);
   ps(:,2:end-1) = tmp(1:4,:);
   xs(:,2:end-1) = tmp(5:7,:);
else
   tmp = reshape(q,[7,length(sol.output_s)]);
   ps = tmp(1:4,:);
   xs = tmp(5:7,:);
    end
end
end
