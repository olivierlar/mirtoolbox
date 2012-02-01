function whilePlaying
global fig
global Fs
global pointerH
global framePos
global player
global CurrentSample
global CurrentFrame
global xlim
global PauseButton
global PPSstate
global smallPointerH
global sliderH
global aH
global followPointerButton

if not(ishandle(fig))
    return
end

CurrentSample=get(player,'CurrentSample');
if CurrentSample<=1
    return
elseif CurrentSample>player.TotalSamples
    set(PauseButton,'Value',1);
    PPSstate='pause';
    stop(player);
    CurrentSample=1;
    CurrentFrame=1;
    %set(pointerH,'XData',xlim(2)+[-.5*frameLength*ones(2,1),zeros(2,1)]); %place pointer in the end of the song (good when followPointer option is activated)
    %drawnow
    return
else
    [tmp, CurrentFrame]=min(framePos(1,:)<(xlim(1)+CurrentSample/Fs));
    %pointerX=framePos([1,2,2,1],CurrentFrame); %[framePos(1,CurrentFrame), framePos(2,CurrentFrame), framePos(2,CurrentFrame), framePos(1,CurrentFrame)];

    set(pointerH,'XData',framePos([1,2,2,1],CurrentFrame));
    set(smallPointerH,'XData',(CurrentSample-1)/player.TotalSamples*[1,1]);
    
    sliderHLim=get(sliderH,'XData');
    sliderWidth=sliderHLim(2)-sliderHLim(1);
    
    if ((CurrentSample-1)/player.TotalSamples> (sliderHLim(1)+sliderWidth*.75) || (CurrentSample-1)/player.TotalSamples< sliderHLim(1)) && get(followPointerButton,'Value') % follow pointer
        
        sliderHLim([1,4])=max(0,min((CurrentSample-1)/player.TotalSamples-sliderWidth*.25,1-sliderWidth));
        sliderHLim([2,3])=min(1,sliderHLim(1)+sliderWidth);
        set(sliderH,'XData',sliderHLim);
        set(aH,'Xlim',xlim(1)+sliderHLim(1:2)*(xlim(2)-xlim(1)));
        
    end
    
    %set(sliderAx)
    
    drawnow
end

end
