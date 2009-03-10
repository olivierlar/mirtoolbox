function val = get(a, propName)
% GET Get properties from the MIRchromagram object
% and return the value

switch propName
    case 'Bins'
        val = get(mirdata(a),'Pos');
    case 'Weight'
        val = get(mirdata(a),'Data');
    otherwise
        val = get(mirdata(a),propName);
end