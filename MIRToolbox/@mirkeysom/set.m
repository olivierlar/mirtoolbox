function cc = set(c,varargin)
% SET Set properties for the MIRchromagram object
% and return the updated object

propertyArgIn = varargin;
d = mirdata(c);
d = set(d,'Title',get(c,'Title'),'Abs',get(c,'Abs'),'Ord',get(c,'Ord'));
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'Weight'
           d = set(d,'Data',val);
       otherwise
           d = set(d,prop,val);
   end
end
cc = class(struct,'mirkeysom',d);