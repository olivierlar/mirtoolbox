function mm = set(m,varargin)
% SET Set properties for the MIRnovelty object
% and return the updated object

propertyArgIn = varargin;
d = mirdata(m);
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   d = set(d,prop,val);
end
mm.method = m.method;
mm = class(mm,'mirnovelty',d);