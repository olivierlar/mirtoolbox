function display(n,ax,songs)
% MIRNOVELTY/DISPLAY display of a novelty curve or diagram

if strcmpi(n.method,'Foote')
    display(mirscalar(n),ax,songs);
    return
end

disp(' ');
d = get(n,'Data');
nam = get(n,'Name');
fp = get(n,'FramePos');
%cha = get(n,'Channels');

if nargin<3 || isempty(songs)
    songs=1:length(d);
end

for i = 1:length(songs)  %For each audio file
    if nargin<2 || isempty(ax)
        figure
    else
        axes(ax)
    end

    fpi = cell2mat(fp{i});
    if size(fpi,1) == 2
        fpi = (fpi(1,:,:,:)+fpi(2,:,:,:))/2;
    end
        
    ng = length(d{i}{1});
    v = zeros(size(fpi,1),ng);
    for j = 1:ng
        v(j,d{i}{1}{j}) = 1;
    end        

    h = imagesc(fpi,1:ng,v);

    title('Novelty graph')
    xlabel('temporal location of frame centers (in s.)')
    ylabel('Granularity levels')
    fig = get(0,'CurrentFigure');
    disp(['The Novelty diagram related to file ',nam{i},' is displayed in Figure ',num2str(fig),'.']);
end
disp(' ');
drawnow