function display(m)
% DISSIMATRIX/DISPLAY display of a dissimilarity matrix
disp(' ');
d = get(m,'Data');
n = get(m,'Name');
fp = get(m,'FramePos');
cha = get(m,'Channels');
t = get(m,'Title');
pp = get(m,'PeakPos');
ap = get(m,'AttackPos');
rp = get(m,'ReleasePos');
for i = 1:length(d)
    if iscell(d{i})
        d{i} = d{i}{1};
    end
    l = size(d{1},3);     % Number of channels
    il = (1-0.15)/l;
    figure
    for k = 1:l         % For each channel
        if l>1
            subplot('Position',[0.1 (k-1)*il+0.1 0.89 il-0.02])
        end
        fpi = cell2mat(fp{i});
        if size(fpi,1) == 2
            fpi = (fpi(1,:,:,:)+fpi(2,:,:,:))/2;
        end
        if strcmp(m.view,'h')
            h = imagesc(fpi,fpi-fpi(end)/2,d{i}(:,:,k));
        else
            h = imagesc(fpi,fpi,d{i}(:,:,k));
        end
        if not(isempty(ap))
            hold on
            for j = 1:length(ap{i}{1})
                if not(isempty(ap{i}{1}{j}))
                    plot([fpi(j);fpi(j)],fpi([ap{i}{1}{j};rp{i}{1}{j}]),'--k')
                    plot(fpi(j),fpi(ap{i}{1}{j}),'dk')
                    plot(fpi(j),fpi(rp{i}{1}{j}),'dk')
                end
            end
        end
        set(gca,'YDir','normal')
        if k == l
            title(t)
        end
        if k == 1
            if strcmp(m.view,'l')
                xlabel('temporal lag (in s.)')
            else
                xlabel('temporal location of frame centers (in s.)')
            end
        end
        if k == ceil(l/2)
            if strcmp(m.view,'h')
                ylabel('relative distance between compared frames (in s.)')
            else
                ylabel('temporal location of frame centers (in s.)')
            end
        end
        if l > 1
            pos = get(gca,'Position');
            hfig = axes('Position',[pos(1)-.05 pos(2)+pos(4)/2 .01 .01],...
                        'Visible','off');
            text(0,0,num2str(cha{i}(k)),'FontSize',12,'Color','r')
        end
    end
    fig = get(0,'CurrentFigure');
    disp(['The ',t,' related to file ',n{i},' is displayed in Figure ',num2str(fig),'.']);
end
disp(' ');
drawnow