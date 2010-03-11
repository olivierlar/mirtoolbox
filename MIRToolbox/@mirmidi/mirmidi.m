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
if not(isamir(x,'mirmidi'))
    o = mironsets(x);
    x = {o x};
end
type = 'mirmidi';
    

function m = main(x,option,postoption)
if iscell(x) %not(isamir(x,'mirmidi'))
    o = x{1};
    a = x{2};
    s = mirsegment(a,o);
    x = mirpitch(s,'Contrast',option.thr);
end
d = get(x,'Data');
fp = get(x,'FramePos');
nmat = cell(1,length(d));
for i = 1:length(d)
    nmat{i} = [];
    for j = 1:length(d{i})
        tij = fp{i}{j}(1);
        dij = diff(fp{i}{j});
        for k = 1:size(d{i}{j},3)
            for l = 1:size(d{i}{j},2)
                for n = 1:length(d{i}{j}{1,l,k})
                    f = d{i}{j}{1,l,k}(n);
                    p = round(hz2midi(f));
                    nmat{i} = [nmat{i}; tij dij 1 p 120 tij dij];
                end
            end
        end
    end
end
m.data = nmat;
m = class(m,'mirmidi',mirdata(x));