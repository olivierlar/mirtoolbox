function cc = set(c,varargin)
% SET Set properties for the MIRchromagram object
% and return the updated object

propertyArgIn = varargin;
pl = c.plabel;
wr = c.wrap;
cl = c.chromaclass;
cf = c.chromafreq;
or = c.register;
d = mirdata(c);
d = set(d,'Title',get(c,'Title'),'Abs',get(c,'Abs'),'Ord',get(c,'Ord'));
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'Magnitude'
           d = set(d,'Data',val);
       case 'Chroma'
           d = set(d,'Pos',val);
       case 'ChromaClass'
           cl = val;
       case 'ChromaFreq'
           cf = val;
       case 'Register'
           or = val;
       case 'PitchLabel'
           pl = val;
       case 'Wrap'
           wr = val;
       otherwise
           d = set(d,prop,val);
   end
end
cc.plabel = pl;
cc.wrap = wr;
cc.chromaclass = cl;
cc.chromafreq = cf;
cc.register = or;
cc = class(cc,'mirchromagram',d);