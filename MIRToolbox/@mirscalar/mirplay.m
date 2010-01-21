function varargout = mirplay(a,varargin)

% mirplay method for mirscalar objects. Help displayed in ../mirplay.m

        ch.key = 'Channel';
        ch.type = 'Integer';
        ch.default = 0;
    option.ch = ch;
        
        sg.key = 'Segment';
        sg.type = 'Integer';
        sg.default = 0;
    option.sg = sg;
        
        se.key = 'Sequence';
        se.type = 'Integer';
        se.default = 0;
    option.se = se;

        inc.key = 'Increasing';
        inc.type = 'MIRtb';
    option.inc = inc;
    
        dec.key = 'Decreasing';
        dec.type = 'MIRtb';
    option.dec = dec;

        every.key = 'Every';
        every.type = 'Integer';
    option.every = every;

specif.option = option;

specif.eachchunk = 'Normal';

varargout = mirfunction(@mirplay,a,varargin,nargout,specif,@init,@main);
if nargout == 0
    varargout = {};
end


function [x type] = init(x,option)
type = '';


function noargout = main(a,option,postoption)
if iscell(a)
    a = a{1};
end
d = get(a,'Data');
if isa(a,'mirpitch')
    amp = get(a,'Amplitude');
end
f = get(a,'Sampling');
n = get(a,'Name');
c = get(a,'Channels');
fp = get(a,'FramePos');
if not(option.se)
    if length(d)>1
        if isfield(option,'inc')
            [unused order] = sort(mirgetdata(option.inc));
        elseif isfield(option,'dec')
            [unused order] = sort(mirgetdata(option.dec),'descend');
        else
            order = 1:length(d);
        end
        if isfield(option,'every')
            order = order(1:option.every:end);
        end
    else
        order = 1;
    end
else
    order = option.se;
end
if not(isempty(order))
    for k = order(:)'
        display(['Playing analysis of file: ' n{k}])   
        dk = d{k};
        if not(iscell(dk))
            dk = {dk};
        end
        if option.ch
            if isempty(c{k})
                chk = option.ch;
            else
                [unused unused chk] = intersect(option.ch,c{k});
            end
        else
            chk = 1:size(dk{1},3);
        end
        if isempty(chk)
            display('No channel to play.');
        end
        for l = chk
            if chk(end)>1
                display(['  Playing channel #' num2str(l)]);
            end
            if option.sg
                sgk = option.sg(find(option.sg<=length(dk)));
            else
                sgk = 1:length(dk);
            end
            for i = sgk
                if sgk(end)>1
                    display(['      Playing segment #' num2str(i)])
                end
                di = dk{i};
                if isa(a,'mirpitch')
                    ampi = amp{k}{i};
                end
                synth = zeros(1,ceil((fp{k}{i}(end)-fp{k}{i}(1))*44100)+1);
                for j = 1:size(di,2)
                    if iscell(di)
                        dj = di{j};
                    else
                        dj = di(:,j);
                    end
                    dj(isnan(dj)) = 0;
                    ampj = zeros(size(dj));
                    if iscell(ampi)
                        ampj(1:size(ampi{j})) = ampi{j};
                    else
                        ampj(1:size(ampi(:,j))) = ampi(:,j);
                    end
                    if not(isempty(dj))
                        k1 = floor((fp{k}{i}(1,j)-fp{k}{i}(1))*44100)+1;
                        k2 = floor((fp{k}{i}(2,j)-fp{k}{i}(1))*44100)+1;
                        if isa(a,'mirpitch')
                            ampj = repmat(ampj,1,k2-k1+1);
                        else
                            ampj = ones(size(dj),k2-k1+1);
                        end
                        synth(k1:k2) = synth(k1:k2) ...
                            + sum(ampj.*sin(2*pi*dj*(0:k2-k1)/44100),1) ...
                                    .*hann(k2-k1+1)';
                        %plot((ampj.*sin(2*pi*dj*(0:k2-k1)/44100))')
                        %drawnow
                    end
                end
                soundsc(synth,44100);
                if sgk(end)>1
                    sound(rand(1,10))
                end
                %pause(0.5)
            end
        end
    end
end
noargout = {};