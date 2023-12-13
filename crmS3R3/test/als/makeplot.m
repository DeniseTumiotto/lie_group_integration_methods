function solcell = makeplot(pattern, refpattern, varargin)

varargin = union(varargin, {});

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
solcell = calc_errors(solcell, refpattern);

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
   figure('Name',[absrel '__' figname]);
   ax = axes();

   if makecpuplot
      xdata = catsolcell(solcell, 'stats.cpu_time');
   else
      % The array of N is needed
      xdata = catsolcell(solcell, 'n');
   end
   
   alsosort = @(x,I) x(I);
   [xdata,I] = sort(xdata);
   
   if hasl
      if hase && absreli == 1
         loglog(ax, ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.q' ]),I), '-o', ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.v' ]),I), '-x', ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.vd']),I), '-v', ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.l' ]),I), '-d', ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.e' ]),I), '-p');
         legend('$q$','$v$','$\dot v$','$\lambda$','$\eta$');
      else
         loglog(ax, ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.q' ]),I), '-o', ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.v' ]),I), '-x', ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.vd']),I), '-v', ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.l' ]),I), '-d');
         legend('$q$','$v$','$\dot v$','$\lambda$');
      end
   else
      loglog(ax, ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.q' ]),I), '-o', ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.v' ]),I), '-x', ...
            xdata, alsosort(catsolcell(solcell, ['err.' absrel '.vd']),I), '-v');
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
   title(strrep([absrel ': ' figname],'_',', '));
end