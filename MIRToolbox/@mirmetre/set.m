function mm = set(m,varargin)
% SET Set properties for the MIRmetre object
% and return the updated object

propertyArgIn = varargin;
ac = m.autocor;
gbpm = m.globpm;
d = mirdata(m);
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'Autocor'
           ac = val;
       case 'Globpm'
           gbpm = val;
       otherwise
           d = set(d,prop,val);
   end
end
mm.autocor = ac;
mm.globpm = gbpm;
mm = class(mm,'mirmetre',d);