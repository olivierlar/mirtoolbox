function val = get(a, propName)
% GET Get properties from the MIRdata object
% and return the value

switch propName
    case 'Mode'
        val = a.mode;
    case 'Legend'
        val = a.legend;
    case 'Parameter'
        val = a.parameter;
    otherwise
        val = get(mirdata(a),propName);
end