function display(m)
% MIRMETRE/DISPLAY display of metrical hierarchy

d = get(m,'Meters');
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
            text(d{h}{1}{i}(i2).timidx(1)-.3,60./d{h}{1}{i}(i2).bpms(1),...
                 num2str(d{h}{1}{i}(i2).lvl))
            for i3 = 1:length(d{h}{1}{i}(i2).score)
                if i3>1
                    plot(d{h}{1}{i}(i2).timidx([i3-1 i3]),...
                        60./d{h}{1}{i}(i2).bpms([i3-1 i3]),'-','Color',irgb);
                end
            end
            for i3 = 1:length(d{h}{1}{i}(i2).score)
                rgb = ones(1,1,3) - ...
                    (d{h}{1}{i}(i2).score(i3) - mic) / micmac * irgb;
                plot(d{h}{1}{i}(i2).timidx(i3),...
                     60./d{h}{1}{i}(i2).bpms(i3),'+','Color',rgb);

            end
        end
    end
end
title('Metrical Hierarchy')