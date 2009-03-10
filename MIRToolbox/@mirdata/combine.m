function c = combine(varargin)

c = varargin{1};
l = length(varargin);
p = cell(1,l);
ch = cell(1,l);
d = cell(1,l);
fp = cell(1,l);
sr = cell(1,l);
n = cell(1,l);
la = cell(1,l);
cl = cell(1,l);
pp = cell(1,l);
pm = cell(1,l);
pv = cell(1,l);
ppp = cell(1,l);
ppv = cell(1,l);
if isa(c,'temporal')
    nb = cell(1,l);
end
if isa(c,'mirscalar')
    m = cell(1,l);
end
for i = 1:l
    argin = varargin{i};
    p{i} = getargin(argin,'Pos');
    ch{i} = getargin(argin,'Channels');
    d{i} = getargin(argin,'Data');
    fp{i} = getargin(argin,'FramePos');
    sr{i} = getargin(argin,'Sampling');
    nb{i} = getargin(argin,'NBits');
    n{i} = getargin(argin,'Name');
    la{i} = getargin(argin,'Label');
    cl{i} = getargin(argin,'Clusters');
    pp{i} = getargin(argin,'PeakPos');
    pm{i} = getargin(argin,'PeakMode');
    pv{i} = getargin(argin,'PeakVal');
    ppp{i} = getargin(argin,'PeakPrecisePos');
    ppv{i} = getargin(argin,'PeakPreciseVal');
    if isa(c,'temporal')
        ct = getargin(argin,'Centered');
        nb{i} = getargin(argin,'NBits');
    end
    if isa(c,'mirscalar')
        m{i} = getargin(argin,'Mode');
    end
end
c = set(c,'Pos',p,'Data',d,'FramePos',fp,'Channels',ch,...
          'Sampling',sr,'NBits',nb,'Name',n,'Label',la,...
          'Clusters',cl,'PeakPos',pp,'PeakMode',pm,'PeakVal',pv,...
          'PeakPrecisePos',ppp,'PeakPreciseVal',ppv);
if isa(c,'temporal')
    c = set(c,'Centered',ct,'NBits',nb);
end
if isa(c,'mirscalar')
    c = set(c,'Mode',m);
end
      
      
function y = getargin(argin,field)
yi = get(argin,field);
if isempty(yi) || ischar(yi)
    y = yi;
else
    y = yi{1};
end