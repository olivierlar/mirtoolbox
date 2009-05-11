function c = combine(varargin)

c = varargin{1};
l = length(varargin);
v = cell(1,l);
vf = cell(1,l);
a = cell(1,l);
af = cell(1,l);
sr = cell(1,l);
n = cell(1,l);
la = cell(1,l);
cl = cell(1,l);
pp = cell(1,l);
pm = cell(1,l);
pv = cell(1,l);
ppp = cell(1,l);
ppv = cell(1,l);
for i = 1:l
    argin = varargin{i};
    v{i} = getargin(argin,'Valence');
    vf{i} = getargin(argin,'ValenceFactors');
    a{i} = getargin(argin,'Activity');
    af{i} = getargin(argin,'ActivityFactors');
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
end
c = set(c,'Valence',v,'ValenceFactors',vf,'Activity',a,'ActivityFactors',af,...
          'Sampling',sr,'NBits',nb,'Name',n,'Label',la,...
          'Clusters',cl,'PeakPos',pp,'PeakMode',pm,'PeakVal',pv,...
          'PeakPrecisePos',ppp,'PeakPreciseVal',ppv);
      
      
function y = getargin(argin,field)
yi = get(argin,field);
if isempty(yi) || ischar(yi)
    y = yi;
else
    y = yi{1};
end