function whilePlaying
global fig
global Fs
global aH
global playheads
global framePos
global player
global CurrentSample
global CurrentFrame
global xlim
global PlayPauseButton
global playIcon
global smallPointerH
global sliderH
global featureAxes
global followPointerButton
global frameSummary
global loopingButton


if not(ishandle(fig))
    return
end
CurrentSample=get(player,'CurrentSample');


%fprintf('Current time: %s\n',num2str(CurrentSample/Fs));
%if CurrentSample<1 || CurrentSample>player.TotalSamples
%set(PlayPauseButton,'Tag','play','cdata',playIcon);
%stop(player);
%CurrentSample=1;
%CurrentFrame(1:end)=1;

if  CurrentSample>1
    
    set(smallPointerH,'XData',(CurrentSample-1)/player.TotalSamples*[1,1]);
    
    CurrentFrame_tmp=max(1,sum(frameSummary<(xlim(1)+CurrentSample/Fs),2));
    
    gothru=find(CurrentFrame_tmp~=CurrentFrame);
    CurrentFrame=CurrentFrame_tmp;
    
    
    for plInd=1:length(gothru)
        x=framePos{gothru(plInd)}(:,CurrentFrame(gothru(plInd)));
        set(playheads{gothru(plInd)},'xdata',getxdata(x));
    end
    
    drawnow
    
    sliderHLim=get(sliderH,'XData');
    sliderWidth=sliderHLim(2)-sliderHLim(1);
    
    if get(loopingButton,'Value')
        
        CurrentSelection=round(sliderHLim([1,2])*player.TotalSamples);
        CurrentSelection(1)=max(1,CurrentSelection(1));
        CurrentSelection(2)=min(player.TotalSamples,CurrentSelection(2));        
        if CurrentSample>=(CurrentSelection(2)-player.TimerPeriod*player.SampleRate)        
            pause(player);
            play(player,CurrentSelection(1));
        end
        
    end
    
    
    if get(followPointerButton,'Value') && ((CurrentSample-1)/player.TotalSamples > (sliderHLim(1)+sliderWidth*.75) || (CurrentSample-1)/player.TotalSamples< sliderHLim(1)) % follow pointer
        
        sliderHLim([1,4])=max(0,min((CurrentSample-1)/player.TotalSamples-sliderWidth*.25,1-sliderWidth));
        sliderHLim([2,3])=min(1,sliderHLim(1)+sliderWidth);
        set(sliderH,'XData',sliderHLim);
        set(aH,'Xlim',xlim(1)+sliderHLim(1:2)*(xlim(2)-xlim(1)));
        for i=1:length(featureAxes)
            set(featureAxes{i},'Xlim',xlim(1)+sliderHLim(1:2)*(xlim(2)-xlim(1)));
        end
        
    end
    
    
end

    function res = getxdata(xdata)
        s=xdata(2)-xdata(1);
        m=mean(xdata);
        res=[xdata([1,2,2])',m+.05*s,m,m-.05*s,xdata([1,1])'];
    end

end
