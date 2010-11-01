function varargout = mirmidi(orig,varargin)
%   m = mirmidi(x) converts into a MIDI sequence.
%   Option associated to mirpitch function can be specified:
%       'Contrast' with default value c = .3

    thr.key = 'Contrast';
    thr.type = 'Integer';
    thr.default = .3;
option.thr = thr;

    mono.key = 'Mono';
    mono.type = 'Boolean';
    mono.default = 1;
option.mono = mono;

specif.option = option;

varargout = mirfunction(@mirmidi,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
try
    hz2midi(440);
catch
    mirerror('MIRMIDI','MIDItoolbox does not seem to be installed.');
end
if not(isamir(x,'mirmidi')) && not(isamir(x,'mirpitch'))
    if isa(x,'mirdesign') && not(option.mono)
        x = set(x,'SeparateChannels',1);
    end
    o = mironsets(x,'Attacks');
    x = {o x};
end
type = 'mirmidi';
    

function m = main(x,option,postoption)
if iscell(x) %not(isamir(x,'mirmidi'))
    o = x{1};
    a = x{2};
    s = mirsegment(a,o);
    x = mirpitch(s,'Contrast',option.thr,'Sum',0);
    do = get(o,'PeakVal');
else
    do = NaN;
end
dp = get(x,'Data');
fp = get(x,'FramePos');
nmat = cell(1,length(dp));
for i = 1:length(dp)
    nmat{i} = [];
    for j = 2:length(dp{i})
        tij = fp{i}{j}(1);
        dij = diff(fp{i}{j});
        if not(iscell(do))
            vij = 120;
        else
            vij = round(do{i}{1}{1}(j-1)/max(do{i}{1}{1})*120);
        end
        for k = 1:size(dp{i}{j},3)
            for l = 1:size(dp{i}{j},2)
                for n = 1:length(dp{i}{j}{1,l,k})
                    f = dp{i}{j}{1,l,k}(n);
                    p = round(hz2midi(f));
                    nmat{i} = [nmat{i}; tij dij 1 p vij tij dij];
                end
            end
        end
    end
end
m = class(struct,'mirmidi',mirdata(x));
m = set(m,'Data',nmat);