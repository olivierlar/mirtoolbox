function ee = set(e,varargin)
% SET Set properties for the MIRemotion object
% and return the updated object

propertyArgIn = varargin;
af = e.activity_fact;
vf = e.valence_fact;
tf = e.tension_fact;
d = mirdata(e);
d = set(d,'Title',get(e,'Title'),'Abs',get(e,'Abs'),'Ord',get(e,'Ord'));
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'ActivityFactors'
           af = val;
       case 'ValenceFactors'
           vf = val;
       case 'TensionFactors'
           tf = val;
       case 'Concepts'
           d = set(d,'Pos',val);
       otherwise
           d = set(d,prop,val);
   end
end
ee.activity_fact = af;
ee.valence_fact = vf;
ee.tension_fact = tf;
ee = class(ee,'miremotion',d);