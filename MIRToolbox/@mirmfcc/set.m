function ss = set(s,varargin)
% SET Set properties for the MIRspectrum object
% and return the updated object

propertyArgIn = varargin;
dl = s.delta;
d = mirdata(s);
d = set(d,'Title',get(s,'Title'),'Abs',get(s,'Abs'),'Ord',get(s,'Ord'));
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'Rank'
           d = set(d,'Pos',val);
       case 'Delta'
           df = val;
       otherwise
           d = set(d,prop,val);
   end
end
ss.delta = dl;
ss = class(ss,'mirmfcc',d);