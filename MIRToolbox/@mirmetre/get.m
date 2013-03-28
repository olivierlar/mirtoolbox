function val = get(a, propName)
% GET Get properties from the MIRmetre object
% and return the value

switch propName
    case 'Autocor'
        val = a.autocor;
    case 'Globpm'
        val = a.globpm;
    otherwise
        val = get(mirdata(a),propName);
end