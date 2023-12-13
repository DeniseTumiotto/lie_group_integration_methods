function makexyyyplot(solcell, pattern, xname, ynames, varargin)


% Only keep sols that match the pattern
sc = filter_solcell(solcell,pattern);

% Create beginning of title name (also used for figure name)
titlename = '';
for k=1:length(ynames)
   titlename = [titlename ynames{k}];
   if not(k==length(ynames))
      titlename = [titlename ',!'];
   end
end
titlename = [titlename '!over!' xname];

% Create pattern name
if ~isempty(pattern)
   patternfields = fieldnames(pattern);
   patternname = [];
   for i=1:length(patternfields)
      patternname = [patternname patternfields{i} '=' ...
                 num2str(pattern.(patternfields{i})) '!!'];
   end
   patternname(end-1:end) = [];
   
   titlename = [titlename '!with!' patternname];
end

% Create figure
figure('Name',strrep(titlename,'!','_'));
ax = axes();
ax.XScale = 'log';
ax.YScale = 'log';

markerlist = 'ox+pd';

% Create plots
hold(ax,'on');
for k=1:length(ynames)
   [xdata,sortind] = sort(catsolcell(sc, xname));
   ydata = catsolcell(sc, ynames{k});
   ydata = ydata(sortind);

   loglog(ax, xdata, ydata, ['-' markerlist(mod(k-1,length(markerlist))+1)], varargin{:});
end

legend(ynames,'interpreter','none');
legend('Location','Best');

% Optical things
title(strrep(strrep(titlename,'!!',', '),'!',' '),'interpreter','none');
xlabel(xname,'interpreter','none');
%ylabel(yname,'interpreter','none');

drawnow;

grid off;
grid on;


% function sc = filter_solcell(solcell,pattern)
% if ~isempty(pattern)
%    sc = cell(1,0);
%    for k=1:length(solcell)
%       if structmatch(solcell{k}, pattern)
%          sc{end+1} = solcell{k};
%       end
%    end
% else
%    sc = solcell;
% end
