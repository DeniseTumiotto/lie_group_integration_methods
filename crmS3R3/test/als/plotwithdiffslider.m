function plotwithdiffslider(mat)

f = figure;
ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
h = plot(ax,mat,'.-');

b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',0, 'min',0, 'max',100);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                'String','100','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String','diff','BackgroundColor',bgcolor);
   
b.Callback = @(es,ed) updateplot(ax,mat,round(es.Value));

function updateplot(ax,mat,k)
plot(ax,mydiff(mat,k),'.-');
title(['diff=' num2str(k)]);
drawnow;
             
function out = mydiff(mat,k)
if k==0
   out = mat;
else
   out = diff(mat,k);
end