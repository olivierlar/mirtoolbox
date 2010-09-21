function varargout = mirmidi(orig,varargin)
%   m = mirmidi(x) converts into a MIDI sequence.
%   Option associated to mirpitch function can be specified:
%       'Contrast' with default value c = .3

    thr.key = 'Contrast';
    thr.type = 'Integer';
    thr.default = .3;
option.thr = thr;

specif.option = option;

varargout = mirfunction(@mirmidi,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
try
    hz2midi(440);
catch
    mirerror('MIRMIDI','MIDItoolbox does not seem to be installed.');
end
if not(isamir(x,'mirmidi')) && not(isamir(x,'mirpitch'))
    o = mironsets(x,'Sum',0);
    x = {o x};
end
type = 'mirmidi';
    

function m = main(x,option,postoption)
if iscell(x) %not(isamir(x,'mirmidi'))
    o = x{1};
    a = x{2};
    da = get(a,'Data');
    nchannels = size(da{1}{1},3);
    nmat = cell(1,length(da));
    for c = 1:nchannels
        da = get(a,'Data');
        do = get(o,'Data');
        for i = 1:length(da)
            for j = 1:length(da{i})
                da{i}{j} = da{i}{j}(:,:,c);
                do{i}{j} = do{i}{j}(:,:,c);
            end
        end
        ac = set(a,'Data',da);
        oc = set(o,'Data',do);
        s = mirsegment(ac,oc);
        x = mirpitch(s,'Contrast',option.thr);
        d = get(x,'Data');
        fp = get(x,'FramePos');
        for i = 1:length(d)
            for j = 1:length(d{i})
                tij = fp{i}{j}(1);
                dij = diff(fp{i}{j});
                for k = 1:size(d{i}{j},3)
                    for l = 1:size(d{i}{j},2)
                        for n = 1:length(d{i}{j}{1,l,k})
                            f = d{i}{j}{1,l,k}(n);
                            p = round(hz2midi(f));
                            nmat{i} = [nmat{i}; tij dij c p 120 tij dij];
                        end
                    end
                end
            end
            if c == nchannels
                nmat{i} = sortrows(nmat{i});
            end
        end
        nmat
    end
end
m.data = nmat;
m = class(m,'mirmidi',mirdata(x));