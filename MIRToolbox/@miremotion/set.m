function ee = set(e,varargin)
% SET Set properties for the MIRemotion object
% and return the updated object

propertyArgIn = varargin;
%dim = e.dim;
dimdata = e.dimdata;
%activity_fact = e.activity_fact;
%valence_fact = e.valence_fact;
%tension_fact = e.tension_fact;
%classes = e.class;
classdata = e.classdata;
%happy_fact = e.happy_fact;
%sad_fact = e.sad_fact;
%tender_fact = e.tender_fact;
%anger_fact = e.anger_fact;
%fear_fact = e.fear_fact;
d = mirdata(e);
d = set(d,'Title',get(e,'Title'),'Abs',get(e,'Abs'),'Ord',get(e,'Ord'));
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'DimData'
           dimdata = val;
       case 'ClassData'
           classdata = val;
       otherwise
           d = set(d,prop,val);
   end
end
ee.dim = e.dim;
ee.dimdata = dimdata;
ee.activity_fact = e.activity_fact;
ee.valence_fact = e.valence_fact;
ee.tension_fact = e.tension_fact;
ee.class = e.class;
ee.classdata = classdata;
ee.happy_fact = e.happy_fact;
ee.sad_fact = e.sad_fact;
ee.tender_fact = e.tender_fact;
ee.anger_fact = e.anger_fact;
ee.fear_fact = e.fear_fact;
ee = class(ee,'miremotion',d);