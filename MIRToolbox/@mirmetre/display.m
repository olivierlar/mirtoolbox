function display(m)
% MIRMETRE/DISPLAY display of metrical hierarchy

d = get(m,'Data');
fp = get(m,'FramePos');
for h = 1:length(d)
    figure,hold on
    for i = 1:length(d{h}{1})
        irgb = shiftdim(1-num2col(i),-1);
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
            timidx = d{h}{1}{i}(i2).timidx;
            text(mean(fp{h}{1}(:,timidx(1)))-.3,...
                 60./d{h}{1}{i}(i2).bpms(1),...
                 num2str(d{h}{1}{i}(i2).lvl))
            for i3 = 1:length(d{h}{1}{i}(i2).score)
                if i3>1 && length(d{h}{1}{i}(i2).globpms) >= i3
                    plot(mean(fp{h}{1}(:,timidx([i3-1 i3]))),...
                         60./d{h}{1}{i}(i2).globpms([i3-1 i3]),...
                         '-+r'); %'Color',irgb);
                end
            end
            for i3 = 1:length(d{h}{1}{i}(i2).score)
                rgb = ones(1,1,3) - ...
                    (d{h}{1}{i}(i2).score(i3) - mic) / micmac * irgb;
                plot(mean(fp{h}{1}(:,timidx(i3))),...
                     60./d{h}{1}{i}(i2).bpms(i3),'+','Color',rgb);

            end
        end
    end
    title('Metrical Hierarchy')
    xlabel('Temporal evolution (in s.)')
    ylabel('Pulsation periods (in s.)')
end
