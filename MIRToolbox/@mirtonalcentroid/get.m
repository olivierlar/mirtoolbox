function val = get(a, propName)
% GET Get properties from the MIRchromagram object
% and return the value

switch propName
    case 'Positions'
        val = get(mirdata(a),'Data');
    case 'Dimensions'
        val = get(mirdata(a),'Pos');
    otherwise
        val = get(mirdata(a),propName);
end