function display(m)
% MATRIX/DISPLAY display of a matrix
disp(' ');
d = get(m,'Data');
n = get(m,'Name');
t = get(m,'Title');
for i = 1:length(d)
    for j = 1:length(d{i})
        figure
        h = imagesc(d{i}{j});
        set(gca,'YDir','normal')
        title(t)
        fig = get(0,'CurrentFigure');
        disp(['The ',t,' related to file ',n{i},' is displayed in Figure ',num2str(fig),'.']);
    end
end
disp(' ');
drawnow