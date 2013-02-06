function pp = set(p,varargin)
% SET Set properties for the MIRpitch object
% and return the updated object

propertyArgIn = varargin;
a = p.amplitude;
ps = p.start;
pe = p.end;
pm = p.mean;
deg = p.degrees;
stb = p.stable;
s = mirscalar(p);
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'Amplitude'
           a = val;
       case 'Start'
           ps = val;
       case 'End'
           pe = val;
       case 'Mean'
           pm = val;
       case 'Degrees'
           deg = val;
       case 'Stable'
           stb = val;
       otherwise
           s = set(s,prop,val);
   end
end
pp.amplitude = a;
pp.start = ps;
pp.end = pe;
pp.mean = pm;
pp.degrees = deg;
pp.stable = stb;
pp = class(pp,'mirpitch',s);