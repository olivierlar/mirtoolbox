function pp = set(p,varargin)
% SET Set properties for the MIRpitch object
% and return the updated object

propertyArgIn = varargin;
a = p.amplitude;
s = mirscalar(p);
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'Amplitude'
           a = val;
       otherwise
           s = set(s,prop,val);
   end
end
pp.amplitude = a;
pp = class(pp,'mirpitch',s);