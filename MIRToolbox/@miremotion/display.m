function display(e)
% EMOTION/DISPLAY display of a 2D emotion
disp(' ');
a = get(e,'Activity');
v = get(e,'Valence');
af = get(e,'ActivityFactors');
vf = get(e,'ValenceFactors');
n = get(e,'Name');
t = get(e,'Title');
for i = 1:length(a)
    seq = NaN;
    if length(a{i})>1
        try
            seq = mirsegment(n{i},0:5:1000);
        end
    end
    [f1 f2] = displayeach(1,v{i},a{i},vf{i},af{i},seq);
    for j = 2:length(a{i})
        displayeach(j,v{i},a{i},vf{i},af{i},seq,f1,f2)
    end
    fig = get(0,'CurrentFigure');
    disp(['The ',t,' related to file ',n{i},' is displayed in Figure ',num2str(fig),'.']);
end
disp(' ');
drawnow


function [f1 f2] = displayeach(j,v,a,vf,af,seq,f1,f2)
if nargin<7
    f1 = figure;
    space = imread('space.png','png');
    image([.4 1.1],[1.7 1.1],space)
    axis([.4 1.1 1.1 1.7])
    set(gca,'YDir','normal')
    hold on 
    xlabel('Valence')
    ylabel('Activity')
    title('Emotion')
    if 1
        f2 = figure;
    else
        f2 = 0;
    end
end

figure(f1)
if j == 1
    plot(v{j},a{j},'o','MarkerSize',12,...
                     'MarkerEdgeColor','k','MarkerFaceColor','k')
else
    plot([v{j-1} v{j}],[a{j-1} a{j}],'o-k')
end

if 1
    figure(f2)
    subplot(1,1,1)

    subplot(2,1,1)
    hold on
    set(gca,'YDir','reverse')
    for k = 1:length(vf{j})
        x1 = 0;% min(0,vf{j}(k));
        y1 = k -.3;
        x2 = x1+abs(vf{j}(k));
        y2 = k +.3;
        pcolor([x1 x1;x2 x2],[y1 y2;y1 y2],[1 1; 1 1])
    end
    set(gca,'YTick',1:7,'YTickLabel',...
        {'Key clarity','Mode','Event density','Repetitiveness',...
         'Tempo','Brightness','Roughness'});
    grid on
    title('Valence factors')

    subplot(2,1,2)
    set(gca,'YDir','reverse')
    hold on
    for k = 1:length(af{j})
        x1 = 0; %min(0,af{j}(k));
        y1 = k -.3;
        x2 = x1+abs(af{j}(k));
        y2 = k +.3;
        pcolor([x1 x1;x2 x2],[y1 y2;y1 y2],[1 1; 1 1])
    end
    set(gca,'YTick',1:7,'YTickLabel',...
        {'Spectral flux','Spectral entropy','(Collapsed version)',...
         'Articulation','Repetitiveness','Fluctuation','Pulse clarity'});
    grid on
    title('Activity factors')
end

drawnow
if isa(seq,'miraudio')
    mirplay(seq,'Segment',j);
end