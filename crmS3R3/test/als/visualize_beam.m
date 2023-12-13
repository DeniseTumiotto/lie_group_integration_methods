function visualize_beam(sol)

make_video = false;

set_view = false;

do_edges = false;

if make_video
   vid = VideoWriter('out.avi');
   open(vid);
end

[~, xs] = q_to_ps_xs(sol.rslt.q(:,1), sol);

figure('outerposition',[0 0 1024 768]);

hold on;
plt = plot3(xs(1,:),xs(2,:),xs(3,:),'-o');
if do_edges
   lin1 = plot3(xs(1,1),xs(2,1),xs(3,1),'-r');
   lin2 = plot3(xs(1,end),xs(2,end),xs(3,end),'-r');
end
axis equal;
grid off;
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
title(['t = ' num2str(sol.rslt.t(:,1))]);

plt.Parent.XLim = [min(min(sol.rslt.q(5:7:end,:)))-1,...
                   max(max(sol.rslt.q(5:7:end,:)))+1];
plt.Parent.YLim = [min(min(sol.rslt.q(6:7:end,:)))-1,...
                   max(max(sol.rslt.q(6:7:end,:)))+1];
plt.Parent.ZLim = [min(min(sol.rslt.q(7:7:end,:)))-1,...
                   max(max(sol.rslt.q(7:7:end,:)))+1];


view1 = [-15, 89];                
%view2 = view1;
view2 = [ 15,   0];
                
%plt.Parent.XLim = [min(xs(1,:))-1, max(xs(1,:))+1];
%plt.Parent.YLim = [min(xs(2,:))-1, max(xs(2,:))+1];
%plt.Parent.ZLim = [min(xs(3,:))-1, max(xs(3,:))+1];

for ii = 1:2:length(sol.rslt.t)
   [~, xs] = q_to_ps_xs(sol.rslt.q(:,ii), sol);
   
   plt.XData = xs(1,:);
   plt.YData = xs(2,:);
   plt.ZData = xs(3,:);
   
   if do_edges
      lin1.XData = [lin1.XData,xs(1,1)];
      lin1.YData = [lin1.YData,xs(2,1)];
      lin1.ZData = [lin1.ZData,xs(3,1)];
      %
      lin2.XData = [lin2.XData,xs(1,end)];
      lin2.YData = [lin2.YData,xs(2,end)];
      lin2.ZData = [lin2.ZData,xs(3,end)];
   end
   
   %plt.Parent.XLim = [min(xs(1,:))-1, max(xs(1,:))+1];
   %plt.Parent.YLim = [min(xs(2,:))-1, max(xs(2,:))+1];
   %plt.Parent.ZLim = [min(xs(3,:))-1, max(xs(3,:))+1];
   
   par = (ii-1)/(length(sol.rslt.t)-1);
   if set_view
      view(view2*par + (1-par)*view1);
   end
   
   title(['t = ' sprintf('%-10.5f',sol.rslt.t(ii))]);
   
   drawnow;
   
   if make_video
      frame = getframe(gcf);
      writeVideo(vid,frame);
   end
end

if make_video
   close(vid);
end

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