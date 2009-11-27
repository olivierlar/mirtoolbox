function val = get(a, propName)
% GET Get properties from the MIRemotion object
% and return the value

switch propName
    case 'Dim'
        val = a.dim;
    case 'Class'
        val = a.class;
    case 'DimData'
        val = a.dimdata;
    case 'ClassData'
        val = a.classdata;
    case 'ActivityFactors'
        val = a.activity_fact;
    case 'ValenceFactors'
        val = a.valence_fact;
    case 'TensionFactors'
        val = a.tension_fact;
    case 'HappyFactors'
        val = a.happy_fact;
    case 'SadFactors'
        val = a.sad_fact;
    case 'TenderFactors'
        val = a.tender_fact;
    case 'AngerFactors'
        val = a.anger_fact;
    case 'FearFactors'
        val = a.fear_fact;
    otherwise
        val = get(mirdata(a),propName);
end