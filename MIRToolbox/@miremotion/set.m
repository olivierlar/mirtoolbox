function ee = set(e,varargin)
% SET Set properties for the MIRemotion object
% and return the updated object

propertyArgIn = varargin;
a = e.activity;
af = e.activity_fact;
v = e.valence;
vf = e.valence_fact;
d = mirdata(e);
d = set(d,'Title',get(e,'Title'));
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'Activity'
           a = val;
       case 'Valence'
           v = val;
       case 'ActivityFactors'
           af = val;
       case 'ValenceFactors'
           vf = val;
       otherwise
           d = set(d,prop,val);
   end
end
ee.valence = v;
ee.valence_fact = vf;
ee.activity = a;
ee.activity_fact = af;
ee = class(ee,'miremotion',d);