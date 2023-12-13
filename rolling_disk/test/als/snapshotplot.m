function snapshotplot(sol,timepoints)

disksteps = 50;
disk = [sin(2*pi*(1:disksteps)/disksteps);
        cos(2*pi*(1:disksteps)/disksteps);
        zeros(1,disksteps)          ];
hold on;
mplt = plot3(sol.rslt.q(5,:),sol.rslt.q(6,:),sol.rslt.q(7,:));
for i=length(sol.rslt.t):-1:1
    pointline(:,i) = applyp(sol.rslt.q(1:4,i),[0,1,0]') + sol.rslt.q(5:7,i);
end
plot3(pointline(1,:), pointline(2,:), pointline(3,:));
    

axis equal;
xlim([min(sol.rslt.q(5,:))-1, max(sol.rslt.q(5,:))+1]);
ylim([min(sol.rslt.q(6,:))-1, max(sol.rslt.q(6,:))+1]);
zlim([min(sol.rslt.q(7,:))-1, max(sol.rslt.q(7,:))+1]);

for i=1:length(timepoints)
   [~,ind] = min(abs(sol.rslt.t - timepoints(i)));
   thisdisk = bsxfun(@plus,applyp(sol.rslt.q(1:4,ind),disk), sol.rslt.q(5:7,ind));

   plot3(thisdisk(1,:),thisdisk(2,:),thisdisk(3,:),'LineWidth',5);
end



function out = applyp(p,v)

out = [p(1)^2*v(1,:)  +  p(2)^2*v(1,:) - (p(3)^2  +  p(4)^2)*v(1,:)  +  p(1)*(-2*p(4)*v(2,:)  +  2*p(3)*v(3,:))  +  2*p(2)*(p(3)*v(2,:)  +  p(4)*v(3,:));
       2*p(2)*p(3)*v(1,:)  +  2*p(1)*p(4)*v(1,:)  +  p(1)^2*v(2,:) - p(2)^2*v(2,:)  +  p(3)^2*v(2,:) - p(4)^2*v(2,:) - 2*p(1)*p(2)*v(3,:)  +  2*p(3)*p(4)*v(3,:);
       p(1)*(-2*p(3)*v(1,:)  +  2*p(2)*v(2,:))  +  2*p(4)*(p(2)*v(1,:)  +  p(3)*v(2,:))  +  p(1)^2*v(3,:) - (p(2)^2  +  p(3)^2 - p(4)^2)*v(3,:)];
