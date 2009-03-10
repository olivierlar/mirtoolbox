function val = get(a, propName)
% GET Get properties from the MIRsimatrix object
% and return the value

switch propName
    case 'DiagWidth'
        val = a.diagwidth;
    otherwise
        val = get(mirdata(a),propName);
end