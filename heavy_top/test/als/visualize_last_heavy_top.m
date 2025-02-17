function sol = visualize_last_heavy_top(sol, alt)
% Visualize the solution of the heavy top calculated with the Fortran
% program heavy_top.

% Calculate refX
refX = qpvkp(qconj(sol.p0(end-4:end)), sol.x0(end-3:end));
% Which Lie group was used
switch sol.liegroup
   case 1 % SO(3)
      name = 'SO(3)';
      getpos = @(q) reshape(q, [3 3])*refX';
   case 2 % S³
      name = 'S³';
      getpos = @(q) qpvkp(q,refX);
   case 3 % SO(3) x R³
      if nargin == 2 && alt == 1
         name = 'SO(3)xR³, rot.';
         getpos = @(q) reshape(q(1:9), [3 3])*refX';
      else
         name = 'SO(3)xR³, pos.';
         getpos = @(q) q(10:12);
      end
   case 4 % S³ x R³
      if nargin == 2 && alt == 1
         name = 'S³xR³, rot.';
         getpos = @(q) qpvkp(q(1:4),refX);
      else
         name = 'S³xR³, pos.';
         getpos = @(q) q(5:7);
      end
   case 5 % SE(3)
      if nargin == 2 && alt == 1
         name = 'SE(3), rot.';
         getpos = @(q) reshape(q(1:9),[3 3])*refX';
      else
         name = 'SE(3), pos.';
         getpos = @(q) q(10:12);
      end
   case 6 % S³ |x R³
      if nargin == 2 && alt == 1
         name = 'S³|xR³, rot.';
         getpos = @(q) qpvkp(q(1:4),refX);
      else
         name = 'S³|xR³, pos.';
         getpos = @(q) q(5:7);
      end
   case 7 % UDQ
      name = 'UDQ, rot.';
      getpos = @(q) qpvkp(q(1:4),refX);
end 
      
      

% Get the solution
tt = sol.rslt.t;
XX = zeros(3,size(sol.rslt.q,2));
for i=1:size(sol.rslt.q,2)
   XX(:,i) = getpos(sol.rslt.q(end-7:end,i));
end

% Save XX for later use
sol.rslt.XX = XX;

% Create figure, axes and first plot
figure();
axe = axes();
plt = plot3(axe, XX(1,1), XX(2,1), XX(3,1),'.-','LineWidth',1.5);
grid on;
axis equal;

% Set limits and title
% Set limits and title
axe.XLim = [min(XX(1,:),[],'all')-0.1, max(XX(1,:),[],'all')+0.1];
axe.YLim = [min(XX(2,:),[],'all')-0.1, max(XX(2,:),[],'all')+0.1];
axe.ZLim = [min(XX(3,:),[],'all')-0.1, max(XX(3,:),[],'all')+0.1];
axe.Title.String = [name ': t=', num2str(tt(1),'%.2f')];

% Loop over intermediate time points
if ~exist('alt','var')
    alt = 1;
end
for i=1:alt:length(tt)
   % Update date of the plot
   plt.XData = XX(1,1:i);
   plt.YData = XX(2,1:i);
   plt.ZData = XX(3,1:i);
   % As well as the title
   axe.Title.String = [name ': t=', num2str(tt(i),'%.2f')];
   %
   pause(0.01);
   % Draw now!
   drawnow;
end

function p = qconj(p)
   p(2:4) = - p(2:4);

function rslt = qpvkp(p,v)
   % first, calculate p*[0\\v], keeping in mind that v(0)=0
   r =  [- p(1+1)*v(1) - p(1+2)*v(2) - p(1+3)*v(3), ...
           p(1+0)*v(1) - p(1+3)*v(2) + p(1+2)*v(3), ...
           p(1+3)*v(1) + p(1+0)*v(2) - p(1+1)*v(3), ...
         - p(1+2)*v(1) + p(1+1)*v(2) + p(1+0)*v(3)];
   % now calculate Im(r*conj(p)), keeping in mind that p*v*conj(p) will
   % always yield a vector and directly applying the conj
   rslt = [r(1+1)*p(1+0) - r(1+0)*p(1+1) + r(1+3)*p(1+2) - r(1+2)*p(1+3), ...
           r(1+2)*p(1+0) - r(1+3)*p(1+1) - r(1+0)*p(1+2) + r(1+1)*p(1+3), ...
           r(1+3)*p(1+0) + r(1+2)*p(1+1) - r(1+1)*p(1+2) - r(1+0)*p(1+3)];