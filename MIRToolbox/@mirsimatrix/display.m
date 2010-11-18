function display(m)
% DISSIMATRIX/DISPLAY display of a dissimilarity matrix
disp(' ');
d = get(m,'Data');
n = get(m,'Name');
fp = get(m,'FramePos');
cha = get(m,'Channels');
t = get(m,'Title');
tp = get(m,'TrackPos');
tv = get(m,'TrackVal');
pp = get(m,'PeakPos');
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
        if strcmp(m.view,'l')
            h = imagesc(fpi,fpi(1:size(d{i},1))-fpi(1),d{i}(:,:,k));
        else
            h = imagesc(fpi,fpi,d{i}(:,:,k));
        end
        if not(isempty(tp)) && not(isempty(tp{i}))
            hold on
            for k = 1:size(tp{i}{1}{1},1)
                prej = 0;
                for j = 1:size(tp{i}{1}{1},2)
                    if tv{i}{1}{1}(k,j)
                        if prej% && not(isempty(tp(k,j)))
                            plot([fpi(prej) fpi(j)],...
                                 [fpi(tp{i}{1}{1}(k,prej)) - fpi(1) ...
                                  fpi(tp{i}{1}{1}(k,j)) - fpi(1)],...
                                 'k','LineWidth',1)
                            plot([fpi(prej) fpi(j)],...
                                 [fpi(tp{i}{1}{1}(k,prej)) - fpi(1) ...
                                  fpi(tp{i}{1}{1}(k,j)) - fpi(1)],...
                                 'w+','MarkerSize',10)
                            plot([fpi(prej) fpi(j)],...
                                 [fpi(tp{i}{1}{1}(k,prej)) - fpi(1) ...
                                  fpi(tp{i}{1}{1}(k,j)) - fpi(1)],...
                                 'kx','MarkerSize',10)
                        end
                        prej = j;
                    end
                end
            end
        elseif not(isempty(pp)) && not(isempty(pp{i}))
            hold on
            for k = 1:length(pp{i}{1})
                for j = 1:length(pp{i}{1}{k})
                    plot(fpi(k),fpi(pp{i}{1}{k}(j)) - fpi(1), ...
                                 'w+','MarkerSize',10)
                end
            end
        end
        set(gca,'YDir','normal')
        if k == l
            title(t)
        end
        if k == 1
            xlabel('temporal location of frame centers (in s.)')
        end
        if k == ceil(l/2)
            if strcmp(m.view,'h')
                ylabel('relative distance between compared frames (in s.)')
            elseif strcmp(m.view,'l')
                ylabel('temporal lag (in s.)')
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