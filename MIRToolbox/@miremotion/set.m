function ee = set(e,varargin)
% SET Set properties for the MIRemotion object
% and return the updated object

propertyArgIn = varargin;
%dim = e.dim;
dimdata = e.dimdata;
activity_fact = e.activity_fact;
valence_fact = e.valence_fact;
tension_fact = e.tension_fact;
%classes = e.class;
classdata = e.classdata;
happy_fact = e.happy_fact;
sad_fact = e.sad_fact;
tender_fact = e.tender_fact;
anger_fact = e.anger_fact;
fear_fact = e.fear_fact;
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
        case 'ActivityFactors'
            activity_fact = val;
        case 'ValenceFactors'
            valence_fact = val;
        case 'TensionFactors'
            tension_fact = val;
        case 'HappyFactors'
            happy_fact = val;
        case 'SadFactors'
            sad_fact = val;
        case 'TenderFactors'
            tender_fact = val;
        case 'AngerFactors'
            anger_fact = val;
        case 'FearFactors'
            fear_fact = val;
       otherwise
           d = set(d,prop,val);
   end
end
ee.dim = e.dim;
ee.dimdata = dimdata;
ee.activity_fact = activity_fact;
ee.valence_fact = valence_fact;
ee.tension_fact = tension_fact;
ee.class = e.class;
ee.classdata = classdata;
ee.happy_fact = happy_fact;
ee.sad_fact = sad_fact;
ee.tender_fact = tender_fact;
ee.anger_fact = anger_fact;
ee.fear_fact = fear_fact;
ee = class(ee,'miremotion',d);