function swapplots
% Read the data from all open plots. Then as much new plots as there were
% unique legend entries in the old plots are created, where there are as
% much legend entries as there were open plots.
% All plots should be usual 2D line plots with legends and titles. They are
% usually obtained by
%   plot(XData1, YData1, fmt <more XData, YData and fmt>);
%   legend(entry1 <more entries>);
%   title(title1);
%
% Imagine the plot lines in a matrix where each column represents a plot,
% each line a legend entry and each entry a line in a plot:
%
%        title1  title2  title3
% leg1   line11  line12  line13
% leg2   line21  line22  line23
%
% where we have three plots with two legend entries each. Applying
% swapplots essentially transposes this matrix and gives back
%
%         leg1    leg2
% title1  line11  line21
% title2  line12  line22
% title3  line13  line23
%
% thus we get two plots with three legend entries each in this case.
%
% The legend entries don't have to appear in each plot.

% Initialize cells
legentries = cell(0);
titles = cell(0);
Xes = cell(0,0);
Yes = cell(0,0);

xlog = false;
ylog = false;

% Read all open plots
gr = groot;
for ifig = 1:length(gr.Children)
   fig = gr.Children(ifig);
   
   if length(fig.Children) ~= 2
      error('swapplots:badNumberOfFigChildren',...
         ['There are ' num2str(length(fig.Children)) ' children of '...
          'the figure, but we need exactly two: one for an axe and '...
          'one for the legend']);
   end
   
   % Get the axe and the leg and the title
   if isa(fig.Children(1),'matlab.graphics.axis.Axes')
      axe = fig.Children(1);
      if isa(fig.Children(2),'matlab.graphics.illustration.Legend')
         leg = fig.Children(2);
      else
         error('swapplots:badFigChildren',...
            'Children of the figure are not a legend and an axis');
      end
   elseif isa(fig.Children(2),'matlab.graphics.axis.Axes')
      axe = fig.Children(2);
      if isa(fig.Children(1),'matlab.graphics.illustration.Legend')
         leg = fig.Children(1);
      else
         error('swapplots:badFigChildren',...
            'Children of the figure are not a legend and an axis');
      end
   else
      error('swapplots:badFigChildren',...
            'Children of the figure are not a legend and an axis');
   end
   if length(leg.String) ~= length(axe.Children)
      error('swapplots:noLegendEntryNotEqualAxeChildren',...
         ['There are not exactly as many legend entries as children'...
          ' of the axe']);
   end
   if isempty(axe.Title.String)
      error('swapplots:emptyTitle',...
         'Title of all plots must not be empty');
   else
      titles{end+1} = axe.Title.String;
   end
   
   % See if we need log scales
   if strcmp(axe.XScale, 'log')
      xlog = true;
   end
   if strcmp(axe.YScale, 'log')
      ylog = true;
   end
   
   % Save labels
   try
      xlab = axe.XLabel.String;
   catch ME
      xlab = '';
   end
   try
      ylab = axe.YLabel.String;
   catch ME
      ylab = '';
   end
   
   Xes{1,end+1} = [];
   Yes{1,end+1} = [];
   
   for i=1:length(leg.String)
      memb = ismember(legentries, axe.Children(i).DisplayName);
      
      if any(memb)
         % Legend entry already known
         legi = find(memb,1);
      else
         legentries{end+1} = axe.Children(i).DisplayName;
         legi = length(legentries);
      end
      
      Xes{legi,end} = axe.Children(i).XData;
      Yes{legi,end} = axe.Children(i).YData;
   end
end

close all;

markers = {'-x', '-o', '-v', '-diamond', '-square'};
while length(titles) > length(markers)
   markers = [markers markers];
end
markers = markers(1:length(titles));

% Create new plots
for i=1:length(legentries)
   figure();
   
   plotargs = reshape([Xes(i,:); Yes(i,:); markers],[1 3*size(Xes,2)]);
   % % Delete all args that are empty
   % plotargs(cellfun(@isempty,plotargs)) = [];
   plot(plotargs{:});
   
   % Set scales
   if xlog
      set(gca, 'XScale', 'log');
   end
   if ylog
      set(gca, 'YScale', 'log');
   end
   
   % Set labels
   xlabel(xlab);
   ylabel(ylab);
   
   % % Set legend
   % legend(titles{not(cellfun(@isempty,Xes(i,:)))});
   legend(titles);
   
   % and title
   title(legentries{i});
   
   % Nice grid
   grid off;
   grid on;
end