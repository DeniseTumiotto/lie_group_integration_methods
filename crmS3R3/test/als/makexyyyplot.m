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
N = length(ynames);
color = linspecer(N+1);
for k=1:N
   [xdata,sortind] = sort(catsolcell(sc, xname));
   ydata = catsolcell(sc, ynames{k});
   ydata = ydata(sortind);

   loglog(ax, xdata, ydata, ['-' markerlist(mod(k-1,length(markerlist))+1)], 'MarkerSize',5,...
       'LineWidth',2,'Color',color(end-(k-1),:),'DisplayName',ynames{k}(end));
end
if nargin > 4
    switch varargin{1}
        case 1
            loglog(ax, xdata, xdata   *((ydata(2))/xdata(2)  ), 'k:', 'LineWidth',2,'DisplayName','1^{st} order');
        case 2
            loglog(ax, xdata, xdata.^2*((ydata(2))/xdata(2)^2), 'k:', 'LineWidth',2,'DisplayName','2^{nd} order');
        case 3
            loglog(ax, xdata, xdata.^3*((ydata(2))/xdata(2)^3), 'k:', 'LineWidth',2,'DisplayName','3^{rd} order');
        case 4
            loglog(ax, xdata, xdata.^4*((ydata(2))/xdata(2)^4), 'k:', 'LineWidth',2,'DisplayName','4^{th} order');
        case 5
            loglog(ax, xdata, xdata.^5*((ydata(2))/xdata(2)^5), 'k:', 'LineWidth',2,'DisplayName','5^{th} order');
        otherwise
            disp('Order not available')
    end
end


ax=gca;
ax.FontSize = 16;

% legend(ynames,'interpreter','none','FontSize',16);
legend('FontSize',20);
legend('Location','Best');

% Optical things
% title(strrep(strrep(titlename,'!!',', '),'!',' '),'interpreter','none','FontSize',20);
if strcmp(ynames{1}(1:7), 'err.abs')
    title('Absolute error vs time step size','FontSize',20);
    ylabel('Absolute Error','interpreter','none','FontSize',20);
elseif strcmp(ynames{1}(1:7), 'err.rel')
    title('Relative error vs time step size','FontSize',20);
    ylabel('Relative Error','interpreter','none','FontSize',20);
else
    disp('Plot not known!')
end
xlabel(xname,'interpreter','none','FontSize',20);

grid off;
grid on;
