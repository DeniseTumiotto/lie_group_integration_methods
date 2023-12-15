function simple_plot(s)

[n, m] = size(s);
leg = cell(1);

for i = 1:max(n, m)

    figure(i)
    hold on
    nTime = length(s{i}.rslt.t);
    for j = 1:floor((nTime-1)/10):nTime
        [~, xs] = q_to_ps_xs(s{i}.rslt.q(:,j), s{i});
        % Plot
        plot3(xs(1,:), xs(2,:), xs(3,:), '-o', 'MarkerSize',3);
        % for the legend
        leg{end+1} = ['t=' num2str(s{i}.rslt.t(j))];
    end
    legend(leg{:});
    view(360,0)
end
% 
% acc = s(1:5:end);
% err = s(2:5:end);
% t_1 = s(3:5:end);
% hnw = s(4:5:end);
% mul = s(5:5:end);
% 
% figure('Name','error')
% plot(t_1, err)
% 
% figure('Name','all time step')
% plot(t_1, hnw)
% 
% figure('Name','accepted time step')
% plot(t_1(find(acc==1)), hnw(find(acc==1)))

end

function [ps, xs] = q_to_ps_xs(q,sol)

ps = zeros(4,sol.n+1);
xs = zeros(3,sol.n+1);


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