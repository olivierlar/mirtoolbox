function stopPlaying

global fig
global player
global CurrentSample
global CurrentFrame
global PlayPauseButton
global playIcon
global followPointerButton


if not(ishandle(fig))
    return
end



CurrentSample=get(player,'CurrentSample');
if CurrentSample>=player.TotalSamples || CurrentSample==1
    set(PlayPauseButton,'Tag','play','cdata',playIcon);
    if ~get(followPointerButton,'Value')
        CurrentSample=1;
        CurrentFrame(1:end)=1;
    end
    
end