function val = get(a, propName)
% GET Get properties from the MIRemotion object
% and return the value

switch propName
    case 'Activity'
        val = a.activity;
    case 'Valence'
        val = a.valence;
    otherwise
        val = get(mirdata(a),propName);
end