function val = get(a, propName)
% GET Get properties from the MIRchromagram object
% and return the value

switch propName
    case 'Strength'
        val = get(mirdata(a),'Data');
    case 'Tonic'
        val = get(mirdata(a),'Pos');
    otherwise
        val = get(mirdata(a),propName);
end