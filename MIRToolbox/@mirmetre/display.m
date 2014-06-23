function display(m)
% MIRMETRE/DISPLAY display of metrical hierarchy

d = get(m,'Data');
ac = get(m,'Autocor');
fp = get(m,'FramePos');
for h = 1:length(d)
    if iscell(ac)
        display(ac{h});
    else
        display(ac);
    end
    for i = 1:length(d{h}{1})
        if i == 1
            irgb = [0 0 0];
        else
            continue
            irgb = num2col(i)/2; %shiftdim(1-num2col(i),-1);
        end
        mac = 0;
        mic = 1;
        for i2 = 1:length(d{h}{1}{i})
            for i3 = 1:length(d{h}{1}{i}(i2).score)
                if d{h}{1}{i}(i2).score(i3) > mac
                    mac = d{h}{1}{i}(i2).score(i3);
                end
                if d{h}{1}{i}(i2).score(i3) < mic
                    mic = d{h}{1}{i}(i2).score(i3);
                end
            end
        end
        micmac = mac-mic;
        if ~micmac
            micmac = 1;
        end
        for i2 = 1:length(d{h}{1}{i})
            if isempty(d{h}{1}{i}(i2).globpms)
                continue
            end
            
            if isempty(d{h}{1}{i}(i2).function)
                %continue
                irgb2 = .33+.66*irgb;
            else
                irgb2 = irgb;
            end
            timidx = d{h}{1}{i}(i2).timidx;
            for i3 = 1:length(d{h}{1}{i}(i2).score)
                if i3>1 && length(d{h}{1}{i}(i2).globpms) >= i3
                    plot(fp{h}{1}(1,timidx([i3-1 i3])),...
                         60./d{h}{1}{i}(i2).globpms([i3-1 i3]),...
                         ':','Color',irgb2);
                end
            end
            for i3 = 1:length(d{h}{1}{i}(i2).globpms)
                scor = (d{h}{1}{i}(i2).score(i3) - mic) / micmac +1e-16;
                %rgb = ones(1,1,3) - scor * irgb;
                plot(fp{h}{1}(1,timidx(i3)),...
                     60./d{h}{1}{i}(i2).bpms(i3),'+','Color',irgb2,'MarkerSize',scor*5);
                plot(fp{h}{1}(1,timidx([i3 i3])),...
                    [60./d{h}{1}{i}(i2).globpms(i3) ...
                     60./d{h}{1}{i}(i2).bpms(i3)],'Color',irgb2);
            end
            if isempty(d{h}{1}{i}(i2).function)
                text(fp{h}{1}(1,timidx(1)),...-1,...
                     60./d{h}{1}{i}(i2).globpms(1),...
                     num2str(d{h}{1}{i}(i2).lvl,3),'FontSize',10,...
                     'Color',[0 0 0],...
                     'BackgroundColor',[1 1 1],'Margin',1)
            else
                text(fp{h}{1}(1,timidx(1)),...-1,...
                     60./d{h}{1}{i}(i2).globpms(1),...
                     num2str(d{h}{1}{i}(i2).lvl,3),'FontSize',15,...
                     'Color',[0 0 0],...
                     'BackgroundColor',[1 1 1],'Margin',1)
            end
        end
    end
    title('Metrical Hierarchy')
    xlabel('Temporal evolution (in s.)')
    ylabel('Pulsation periods (in s.)')
end
