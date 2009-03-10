function d = mirdecay(x,thr);
%   d = mirdecay(x) evaluates the decay time (in s). More precisely, returns
%       the temporal distance between the last significant local maximum and the end of the audio. A local maximum is
%       considered as significant if the ratio between the corresponding
%       amplitude at that point and the amplitude of the global maximum
%       exceeds a given threshold thr.
%   Optional argument:
%       d = mirdecay(x,thr) species the value of that significancy threshold.
%           Default value: 0.5

if nargin<2
    thr = 0.5;
end
e = mirenvelope(x);
ve = get(e,'Data');
de = mirenvelope(e,'Diff');
t = get(de,'Time');
vd = get(de,'Data');
dt = cell(1,length(vd));
for h = 1:length(vd)
    dt{h} = cell(1,length(vd{h}));
    for i = 1:length(vd{h})
        ti = t{h}{i};
        vdi = vd{h}{i};
        vei = ve{h}{i};
        z = find(and(vdi(1:end-1)>0,vdi(2:end)<0));
        sz = find(vei(z)/max(vei) > thr);
        dt{h}{i} = ti(end) - ti(z(sz(end)));
    end
end
d = mirscalar(e,'Data',dt,'Title','Decay time','Unit','s.');