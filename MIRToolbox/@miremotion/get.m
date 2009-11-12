function val = get(a, propName)
% GET Get properties from the MIRemotion object
% and return the value

switch propName
    case 'Concepts'
        val = get(mirdata(a),'Pos');
    case 'ActivityFactors'
        val = a.activity_fact;
    case 'ValenceFactors'
        val = a.valence_fact;
    case 'TensionFactors'
        val = a.tension_fact;
    otherwise
        val = get(mirdata(a),propName);
end