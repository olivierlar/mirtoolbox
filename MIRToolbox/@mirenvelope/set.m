function ee = set(e,varargin)
% SET Set properties for the MIRenvelope object
% and return the updated object

propertyArgIn = varargin;
ds = e.downsampl;
hw = e.hwr;
lg = e.log;
df = e.diff;
mt = e.method;
ph = e.phase;
t = mirtemporal(e);
t = set(t,'Title',get(e,'Title'),'Abs',get(e,'Abs'),'Ord',get(e,'Ord'));
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'DownSampling'
           ds = val;
       case 'Halfwave'
           hw = val;
       case 'Diff'
           df = val;
       case 'Log'
           lg = val;
   otherwise
           t = set(t,prop,val);
   end
end
ee.downsampl = ds;
ee.hwr = hw;
ee.diff = df;
ee.log = lg;
ee.method = mt;
ee.phase = ph;
ee = class(ee,'mirenvelope',t);