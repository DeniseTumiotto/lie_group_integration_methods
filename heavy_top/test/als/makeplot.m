function solcell = makeplot(pattern, refpattern, varargin)

varargin = union(varargin, {});

if ismember('SO3xR3', varargin)
   SO3xR3 = true;
   errname = 'SO3xR3_err';
else
   SO3xR3 = false;
   errname = 'err';
end
if ismember('rel', varargin)
   makerelplot = true;
else
   makerelplot = false;
end
if ismember('abs', varargin)
   makeabsplot = true;
else
   makeabsplot = false;
end
if ismember('cpu', varargin)
   makecpuplot = true;
else
   makecpuplot = false;
end

% Load results that match pattern and calculate errors
solcell = load_all_config_and_bin(pattern);
if SO3xR3
   refpattern.liegroup = 3;
   refcell = load_all_config_and_bin(refpattern);
   solcell = calc_SO3xR3_errors(solcell, refcell);
else
   solcell = calc_errors(solcell, refpattern);
end

% Get field names of the pattern
patternfields = fieldnames(pattern);

% Create figure name (also used for title) from the pattern
figname = [];
for i=1:length(patternfields)
   figname = [figname patternfields{i} '=' ...
              num2str(pattern.(patternfields{i})) '_'];
end
figname(end) = [];

% Find out if there is some lamba
hasl = false;
for i=1:length(solcell)
   if isfield(solcell{i}.rslt,'l')
      hasl = true;
   end
end

% Find out if there is some eta
hase = false;
for i=1:length(solcell)
   if isfield(solcell{i}.rslt,'e')
      hase = true;
   end
end

% Cell that holds abs and rel
absrelcell = {};
if makeabsplot
   absrelcell{end+1} = 'abs';
end
if makerelplot
   absrelcell{end+1} = 'rel';
end

for absreli = 1:length(absrelcell)
   absrel = absrelcell{absreli};
   
   % Create figure
   figure('Name',[absrel '_' errname '__' figname]);
   ax = axes();

   if makecpuplot
      xdata = catsolcell(solcell, 'stats.cpu_time');
   else
      % The array of N is needed
      xdata = catsolcell(solcell, 'steps');
   end
   
   if hasl
      if hase && absreli == 1
         loglog(ax, ...
            xdata, catsolcell(solcell, [errname '.' absrel '.q' ]), '-o', ...
            xdata, catsolcell(solcell, [errname '.' absrel '.v' ]), '-x', ...
            xdata, catsolcell(solcell, [errname '.' absrel '.vd']), '-v', ...
            xdata, catsolcell(solcell, [errname '.' absrel '.l' ]), '-d', ...
            xdata, catsolcell(solcell, [errname '.' absrel '.e' ]), '-p');
         legend('$q$','$v$','$\dot v$','$\lambda$','$\eta$');
      else
         loglog(ax, ...
            xdata, catsolcell(solcell, [errname '.' absrel '.q' ]), '-o', ...
            xdata, catsolcell(solcell, [errname '.' absrel '.v' ]), '-x', ...
            xdata, catsolcell(solcell, [errname '.' absrel '.vd']), '-v', ...
            xdata, catsolcell(solcell, [errname '.' absrel '.l' ]), '-d');
         legend('$q$','$v$','$\dot v$','$\lambda$');
      end
   else
      loglog(ax, ...
            xdata, catsolcell(solcell, [errname '.' absrel '.q' ]), '-o', ...
            xdata, catsolcell(solcell, [errname '.' absrel '.v' ]), '-x', ...
            xdata, catsolcell(solcell, [errname '.' absrel '.vd']), '-v');
      legend('$q$','$v$','$\dot v$');
   end
   
   % Make nicely spaced grid
   grid off;
   grid on;
   
   % Set labels and title
   if makecpuplot
      xlabel('CPU time');
   else
      xlabel('steps');
   end
   ylabel('$maxerr$');
   title(strrep([absrel ' ' errname ': ' figname],'_',', '));
end