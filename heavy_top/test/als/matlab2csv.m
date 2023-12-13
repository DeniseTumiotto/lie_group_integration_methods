function matlab2csv(path, ax, nsteps)

if nargin < 3
   nsteps = 1;
end

if nargin < 2
   ax = gca;
end

if nargin < 1
   path = '.';
end

if not(path(end) == '/')
   path = [path '/'];
end


if ~exist(path, 'dir')
   mkdir(path)
end

if not(strcmp(ax.Type,'axes'))
   error('Object is not an axe');
end

if isempty(ax.Title.String)
   tit = ['fig' num2str(double(ax))];
else
   tit = fnamify(ax.Title.String);
end

% Assuming that the data is a line or scatterplot
line_obj_handles = findobj(ax,'type','line');

for k=1:length(line_obj_handles)
   li = line_obj_handles(k);
   if isempty(li.DisplayName)
      nam = num2str(k);
   else
      nam = fnamify(li.DisplayName);
   end
   % Write data
   csvwrite([path tit '--' nam '.csv'], [li.XData([1:nsteps:end-1 end])', li.YData([1:nsteps:end-1 end])']);
end

function str = fnamify(str)

for k=1:length(str)
   if not((97<=str(k) && str(k)<=122) || ...
          (65<=str(k) && str(k)<=90)  || ...
          (48<=str(k) && str(k)<=57)    )
      str(k) = '_';
   end
end
