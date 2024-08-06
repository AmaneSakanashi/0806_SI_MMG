function SaveFig(fig, filename)
%
% Save figure (arg:fig) to filename.png and .pdf
% do not put extention(.png or .pdf) to arg filename 
%
% example: SaveFig(figure(1),'example')

% draw the figure
drawnow;

% change paper size 
temp.figunit = fig.Units;
fig.Units = 'centimeters';
pos = fig.Position;
fig.PaperPositionMode = 'Auto';
temp.figpaperunit = fig.PaperUnits;
fig.PaperUnits = 'centimeters';
temp.figsize = fig.PaperSize;
fig.PaperSize = [pos(3), pos(4)];
set(gca, 'LooseInset', get(gca, 'TightInset'));
% save to file
filename_pdf = strcat(filename, '.pdf');
filename_png = strcat(filename, '.png');
print(fig,filename_pdf,'-dpdf','-r300','-bestfit')
print(fig,filename_png,'-dpng','-r300')
% resume paper size
fig.PaperSize = temp.figsize;
fig.Units = temp.figunit;
fig.PaperUnits = temp.figpaperunit;

end