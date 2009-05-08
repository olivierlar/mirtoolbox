function display(e)
% EMOTION/DISPLAY display of a 2D emotion
disp(' ');
a = get(e,'Activity');
v = get(e,'Valence');
n = get(e,'Name');
t = get(e,'Title');
space = imread('space.png','png');
for i = 1:length(a)
    for j = 1:length(a{i})
        figure
        image([0 1],[1 0],space)
        set(gca,'YDir','normal')
        hold on    
        %nframes = length(a{i}{j});
        %for k = 1:length(a{i}{j})
        %    siz = 1-(k-1)/nframes;
        %    plot(a{i}{j}(k),v{i}{j}(k),'o','MarkerSize',12*(siz).^(1/3),...
        %            'MarkerEdgeColor','k','MarkerFaceColor',[siz siz siz])
        %end
        plot(v{i}{j}(1),a{i}{j}(1),'o','MarkerSize',12,...
                        'MarkerEdgeColor','k','MarkerFaceColor','k')
        plot(v{i}{j},a{i}{j},'o-k')
        plot(v{i}{j}(end),a{i}{j}(end),'o','MarkerSize',12,...
                        'MarkerEdgeColor','k','MarkerFaceColor','w')
        xlabel('Valence')
        ylabel('Activity')
        title(t)
        fig = get(0,'CurrentFigure');
        disp(['The ',t,' related to file ',n{i},' is displayed in Figure ',num2str(fig),'.']);
    end
end
disp(' ');
drawnow