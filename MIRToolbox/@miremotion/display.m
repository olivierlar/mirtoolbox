function display(e)
% MIREMOTION/DISPLAY display of emotion contents
d = get(e,'Dim');
dd = get(e,'DimData');
c = get(e,'Class');
cd = get(e,'ClassData');
fp = get(e,'FramePos');
n = get(e,'Name');
if 0 %not(isempty(d))
    if length(d)==1
        s = mirscalar(e);
        s = set(s,'Pos',{{d'}},'Data',dd,...
            'Title','Dimensional Emotion','Abs','Emotion Dimensions')
    else
        for i = 1:length(dd)  % For each audio file
            figure
            hold on
            ddi = dd{i};
            fpi = fp{i};
            ddo = NaN;
            xlabel(d{1});
            ylabel(d{2});
            title('Dimensional Emotions');
            t = 0;
            if length(d)==3
                grid on
                axis square
                rotate3d
                view(-46,19)
                zlabel(d{3});
            end
            for j = 1:length(ddi)
                ddj = ddi{j};
                fpj = fpi{j};
                if isnan(ddo)
                    if length(d)==2
                        plot(ddj(1),ddj(2),'+');
                    else
                        plot3(ddj(1),ddj(2),ddj(3),'+');
                    end
                else
                    if length(d)==2
                        plot([ddo(1) ddj(1)],[ddo(2) ddj(2)],'+-');
                    else
                        plot3([ddo(1) ddj(1)],[ddo(2) ddj(2)],...
                             [ddo(3) ddj(3)],'+-');
                    end
                end
                if length(ddi)>1 && mean(fpj)>=t
                    if length(d)==2
                        text(ddj(1),ddj(2),num2str(mean(fpj),'%2.1f'),...
                            'FontSize',15,'HorizontalAlignment','Center');
                    else
                        text(ddj(1),ddj(2),ddj(3),num2str(mean(fpj),'%2.1f'),...
                            'FontSize',15,'HorizontalAlignment','Center');
                    end
                    t = ceil(mean(fpj));
                end
                ddo = ddj;
            end
        end
        display(['The Dimensional Emotion related to file ',n{i},...
            ' is displayed in Figure ',gcf'.']);
    end
end
if not(isempty(d))
    s = mirscalar(e);
    s = set(s,'Pos',{{d'}},'Data',dd,...
        'Title','Dimensional Emotion','Abs','Emotion Dimensions')
end
if not(isempty(c))
    s = mirscalar(e);
    s = set(s,'Pos',{{c'}},'Data',cd,...
        'Title','Basic Emotion Set','Abs','Emotion Classes')
end