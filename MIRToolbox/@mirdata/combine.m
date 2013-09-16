function c = combine(varargin)

c = varargin{1};
l = length(varargin);
p = cell(1,l);
ch = cell(1,l);
d = cell(1,l);
fp = cell(1,l);
fr = cell(1,l);
sr = cell(1,l);
n = cell(1,l);
la = cell(1,l);
le = cell(1,l);
cl = cell(1,l);
pp = cell(1,l);
pm = cell(1,l);
pv = cell(1,l);
ppp = cell(1,l);
ppv = cell(1,l);
tp = cell(1,l);
tv = cell(1,l);
tpp = cell(1,l);
tpv = cell(1,l);
ap = cell(1,l);
rp = cell(1,l);
if isa(c,'temporal')
    nb = cell(1,l);
end
if isa(c,'mirscalar')
    m = cell(1,l);
    if isa(c,'mirpitch')
        pa = cell(1,l);
        ps = cell(1,l);
        pe = cell(1,l);
        pi = cell(1,l);
        pd = cell(1,l);
    end
end
if isa(c,'miremotion')
    dd = cell(1,l);
    cd = cell(1,l);
end
if isa(c,'mirmetre')
    ac = cell(1,l);
    g = cell(1,l);
end
for i = 1:l
    argin = varargin{i};
    p{i} = getargin(argin,'Pos');
    ch{i} = getargin(argin,'Channels');
    d{i} = getargin(argin,'Data');
    fp{i} = getargin(argin,'FramePos');
    fr{i} = getargin(argin,'FrameRate');
    sr{i} = getargin(argin,'Sampling');
    nb{i} = getargin(argin,'NBits');
    n{i} = getargin(argin,'Name');
    la{i} = getargin(argin,'Label');
    le{i} = getargin(argin,'Length');
    cl{i} = getargin(argin,'Clusters');
    pp{i} = getargin(argin,'PeakPos');
    pm{i} = getargin(argin,'PeakMode');
    pv{i} = getargin(argin,'PeakVal');
    ppp{i} = getargin(argin,'PeakPrecisePos');
    ppv{i} = getargin(argin,'PeakPreciseVal');
    tp{i} = getargin(argin,'TrackPos');
    tv{i} = getargin(argin,'TrackVal');
    tpp{i} = getargin(argin,'TrackPrecisePos');
    tpv{i} = getargin(argin,'TrackPreciseVal');
    ap{i} = getargin(argin,'AttackPos');
    rp{i} = getargin(argin,'ReleasePos');
    if isa(c,'temporal')
        ct = getargin(argin,'Centered');
        nb{i} = getargin(argin,'NBits');
    end
    if isa(c,'mirscalar')
        m{i} = getargin(argin,'Mode');
        if isa(c,'mirpitch')
            pa{i} = getargin(argin,'Amplitude');
            ps{i} = getargin(argin,'Start');
            pe{i} = getargin(argin,'End');
            pi{i} = getargin(argin,'Mean');
            pd{i} = getargin(argin,'Degrees');
        end
    end
    if isa(c,'miremotion')
        dd{i} = getargin(argin,'DimData');
        cd{i} = getargin(argin,'ClassData');
    end
    if isa(c,'mirmetre')
        ac{i} = getargin(argin,'Autocor');
        g{i} = getargin(argin,'Globpm');
    end
end
c = set(c,'Pos',p,'Data',d,'FramePos',fp,'FrameRate',fr,'Channels',ch,...
          'Sampling',sr,'NBits',nb,'Name',n,'Label',la,'Length',le,...
          'Clusters',cl,'PeakPos',pp,'PeakMode',pm,'PeakVal',pv,...
          'PeakPrecisePos',ppp,'PeakPreciseVal',ppv,...
          'TrackPos',tp,'TrackVal',tv,...
          'TrackPrecisePos',tpp,'TrackPreciseVal',tpv,...
          'AttackPos',ap,'ReleasePos',rp);
if isa(c,'temporal')
    c = set(c,'Centered',ct,'NBits',nb);
end
if isa(c,'mirscalar')
    c = set(c,'Mode',m);
    if isa(c,'mirpitch')
        c = set(c,'Amplitude',pa,'Start',ps,'End',pe,'Mean',pi,...
                  'Degrees',pd);
    end
end
if isa(c,'miremotion')
    c = set(c,'DimData',dd,'ClassData',cd);
end
if isa(c,'mirmetre')
    c = set(c,'Autocor',ac,'Globpm',g);
end
      
      
function y = getargin(argin,field)
yi = get(argin,field);
if isempty(yi) || ischar(yi) || ~iscell(yi)
    y = yi;
else
    y = yi{1};
end