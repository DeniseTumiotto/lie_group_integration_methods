function makexyplot(solcell, pattern, xname, yname, byname, varargin)


% Only keep sols that match the pattern
sc = filter_solcell(solcell,pattern);

% Create beginning of title name (also used for figure name)
titlename = [yname '!over!' xname];

% Get number of lines and their legend data
if ~isempty(byname)
   bydata = unique(catsolcell(sc, byname));
   titlename = [titlename '!by!' byname];
end

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
if ~isempty(byname)
   hold(ax,'on');
   for k=1:length(bydata)
      bysc = filter_solcell(sc, struct(byname,bydata(k)));
      [xdata,sortind] = sort(catsolcell(bysc, xname));
      ydata = catsolcell(bysc, yname);
      ydata = ydata(sortind);

      loglog(ax, xdata, ydata, ['-' markerlist(mod(k-1,length(markerlist))+1)], varargin{:});

      legendstr{k} = [byname '=' num2str(bydata(k))];
   end

   legend(legendstr,'interpreter','none');
   legend('Location','Best');
else
   xdata = catsolcell(sc, xname);
   ydata = catsolcell(sc, yname);

   plot(ax, xdata, ydata, varargin{:});
end

% Optical things
title(strrep(strrep(titlename,'!!',', '),'!',' '),'interpreter','none');
xlabel(xname,'interpreter','none');
ylabel(yname,'interpreter','none');

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
