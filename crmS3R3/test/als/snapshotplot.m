function snapshotplot(sol, varargin)
figure()
hold on;

if nargin < 1
    error('Not enough arguments!')
else
    for j = 2:2:nargin
        switch lower(varargin{j-1})
            case {'extreme','extr traj'}
                traj_extr = varargin{j};
            case {'stride','number'}
                stride = varargin{j};            
            case {'color','bw'}
                bw = varargin{j};
            otherwise
                warning(['parameter ''' varargin{j-1} ''' not recognized']);
        end
    end
end

if ~exist('traj_extr','var')
    traj_extr = 0;
end
if ~exist('stride','var')
    stride = 8;
end
if ~exist('bw','var')
    bw = 0;
end

ts =  [sol.rslt.t(1:ceil(length(sol.rslt.t)/stride):end) sol.rslt.t(end)];

if traj_extr
    if sol.fixed_x0 == 0 && sol.fixed_p0 == 0
       plot3(sol.rslt.q(5,1:stride:end), sol.rslt.q(6,1:stride:end), sol.rslt.q(7,1:stride:end),...
           ':k', 'LineWidth',2,'DisplayName','Trajectory left end');
    end
    if sol.fixed_xn == 0 && sol.fixed_pn == 0
       plot3(sol.rslt.q(end-2,1:stride:end), sol.rslt.q(end-1,1:stride:end), sol.rslt.q(end,1:stride:end),...
           ':k', 'LineWidth',2,'DisplayName','Trajectory right end');
    end
end

if bw
    color = linspecer(stride+3,'gray');
else
    color = linspecer(stride+1);
end
for it=1:stride+1
   t = ts(it);

   % Minimal index
   [~, ind] = min(abs(sol.rslt.t - t));
   [~, xs] = q_to_ps_xs(sol.rslt.q(:,ind), sol);

   % Plot
   plot3(xs(1,:), xs(2,:), xs(3,:), '-o', 'MarkerSize',5, 'LineWidth',3, ...
         'DisplayName', sprintf('t=%.1f s', t),'MarkerFaceColor', color(end-(it-1),:),...
         'Color', color(end-(it-1),:));
end
ax=gca;
ax.FontSize = 16;
ax.DataAspectRatio = [1,1,1];
box on
if strcmp(sol.problem_name,'roll-up')
    view([0 -1 0])
    ax.XLim = [-4, 10];
    ax.ZLim = [-6, 2];
elseif strcmp(sol.problem_name,'flying_spaghetti')
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
