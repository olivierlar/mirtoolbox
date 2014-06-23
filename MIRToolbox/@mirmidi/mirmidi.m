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

    release.key = {'Release','Releases'};
    release.type = 'String';
    release.choice = {'Olivier','Valeri',0,'no','off'};
    release.default = 'Valeri';
option.release = release;

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
    o = mironsets(x,'Attacks','Releases',option.release);
    x = {o x};
end
type = 'mirmidi';
    

function m = main(x,option,postoption)
transcript = 0;
if iscell(x)
    o = x{1};
    do = get(o,'PeakVal');
    da = get(o,'AttackPosUnit');
    dr = get(o,'ReleasePosUnit');
    a = x{2};
    s = mirsegment(a,o);
    x = mirpitch(s,'Contrast',option.thr,'Sum',0);
    % x = mircentroid(s);
    dp = get(x,'Data');
else
    do = [];
    if isa(x,'mirpitch')
        da = get(x,'Start');
        dr = get(x,'End');
        dp = get(x,'Degrees');
        if isempty(da)
            dp = get(x,'Data');
        else
            transcript = 1;
        end
    else
        da = get(x,'AttackPosUnit');
        dr = get(x,'ReleasePosUnit');
    end
end
df = get(x,'FramePos');
nmat = cell(1,length(dp));
if transcript
    for i = 1:length(dp)
        nmat{i} = zeros(length(dp{i}{1}{1}),7);
        for j = 1:length(dp{i}{1}{1})
            t = df{i}{1}(1,da{i}{1}{1}(j));
            d = df{i}{1}(2,dr{i}{1}{1}(j) )- t;
            v = 120;
            p = dp{i}{1}{1}(j) + 62;
            nmat{i}(j,:) = [t d 1 p v t d];
        end
    end
else
    for i = 1:length(dp)
        nmat{i} = [];
        if isempty(do)
            first = 1;
        else
            first = 2;
        end
        for j = first:length(dp{i})
            if isempty(do)
                tij = df{i}{j}(1);
                dij = df{i}{j}(2)- tij;
                vij = 120;
            else
                tij = da{i}{1}{1}(j-1);
                if isempty(dr{i})
                    dij = 0;
                else
                    dij = dr{i}{1}{1}(j-1) - tij;
                end
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
end
m = class(struct,'mirmidi',mirdata(x));
m = set(m,'Data',nmat);