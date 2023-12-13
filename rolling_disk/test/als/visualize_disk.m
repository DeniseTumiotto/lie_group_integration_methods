function visualize_disk(sol)

do_film = false;

disksteps = 50;
disk = [sin(2*pi*(1:disksteps)/disksteps);
        cos(2*pi*(1:disksteps)/disksteps);
        zeros(1,disksteps)          ];
     
dplt = plot3(disk(1,:),disk(2,:),disk(3,:),'LineWidth',5);
hold on;
mplt = plot3(0,0,1);
hold off;
axis equal;
xlim([min(sol.rslt.q(5,:))-1, max(sol.rslt.q(5,:))+1]);
ylim([min(sol.rslt.q(6,:))-1, max(sol.rslt.q(6,:))+1]);
zlim([min(sol.rslt.q(7,:))-1, max(sol.rslt.q(7,:))+1]);
%view([-90,0])
     
if do_film
   v = VideoWriter('output.avi');
   open(v);
   
   fig = gcf;
   fig.Position(3) = 1024;
   fig.Position(4) = 768;
end

for i=1:1:length(sol.rslt.t)
   thisdisk = bsxfun(@plus,applyp(sol.rslt.q(1:4,i),disk), sol.rslt.q(5:7,i));
   dplt.XData = thisdisk(1,:);
   dplt.YData = thisdisk(2,:);
   dplt.ZData = thisdisk(3,:);
   title(['t = ' num2str(sol.rslt.t(i),'%.3f')])
   
   set(mplt,'XData', [mplt.XData sol.rslt.q(5,i)],...
            'YData', [mplt.YData sol.rslt.q(6,i)],...
            'ZData', [mplt.ZData sol.rslt.q(7,i)]);
   
   drawnow;
   
   if do_film
      frame = getframe(gcf);
      writeVideo(v,frame);
      if i==1 || i==length(sol.rslt.t)
         for k=1:15
         writeVideo(v,frame);
         end
      end
   end
   pause(0.02);
   
end

if do_film
   close(v);
end
     
     
     
function out = applyp(p,v)

out = [p(1)^2*v(1,:)  +  p(2)^2*v(1,:) - (p(3)^2  +  p(4)^2)*v(1,:)  +  p(1)*(-2*p(4)*v(2,:)  +  2*p(3)*v(3,:))  +  2*p(2)*(p(3)*v(2,:)  +  p(4)*v(3,:));
       2*p(2)*p(3)*v(1,:)  +  2*p(1)*p(4)*v(1,:)  +  p(1)^2*v(2,:) - p(2)^2*v(2,:)  +  p(3)^2*v(2,:) - p(4)^2*v(2,:) - 2*p(1)*p(2)*v(3,:)  +  2*p(3)*p(4)*v(3,:);
       p(1)*(-2*p(3)*v(1,:)  +  2*p(2)*v(2,:))  +  2*p(4)*(p(2)*v(1,:)  +  p(3)*v(2,:))  +  p(1)^2*v(3,:) - (p(2)^2  +  p(3)^2 - p(4)^2)*v(3,:)];