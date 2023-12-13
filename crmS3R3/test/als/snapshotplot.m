function snapshotplot(sol, ts, stride)

hold on;

if nargin < 2
   %ts =  0:2.5:15;
   ts = [0 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
end

if nargin < 3
   stride = 8;
end

leg = {};

if sol.fixed_x0 == 0 && sol.fixed_p0 == 0
   plot3(sol.rslt.q(5,1:stride:end), sol.rslt.q(6,1:stride:end), sol.rslt.q(7,1:stride:end),'b');
   leg{end+1} = 'trajectory_left_end';
end
if sol.fixed_xn == 0 && sol.fixed_pn == 0
   plot3(sol.rslt.q(end-2,1:stride:end), sol.rslt.q(end-1,1:stride:end), sol.rslt.q(end,1:stride:end),'b');
   leg{end+1} = 'trajectory_right_end';
end



for it=1:length(ts)
   t = ts(it);

   % Minimal index
   [~, ind] = min(abs(sol.rslt.t - t));

   [~, xs] = q_to_ps_xs(sol.rslt.q(:,ind), sol);

   % Plot
   plot3(xs(1,:), xs(2,:), xs(3,:), 'r-o', 'MarkerSize',3);

   % for the legend
   leg{end+1} = ['t=' num2str(t)];
end

legend(leg{:});


xlabel('x');
ylabel('y');
zlabel('z');

grid off;
grid on;

axis equal;
axis tight;

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
