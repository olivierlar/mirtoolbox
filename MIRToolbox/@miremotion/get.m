function val = get(a, propName)
% GET Get properties from the MIRemotion object
% and return the value

switch propName
    case 'Activity'
        val = a.activity;
    case 'Valence'
        val = a.valence;
    case 'ActivityFactors'
        val = a.activity_fact;
    case 'ValenceFactors'
        val = a.valence_fact;
    otherwise
        val = get(mirdata(a),propName);
end