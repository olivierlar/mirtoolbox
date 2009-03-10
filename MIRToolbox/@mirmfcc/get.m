function val = get(a, propName)
% GET Get properties from the MIRspectrum object
% and return the value

switch propName
    case 'Rank'
        val = get(mirdata(a),'Pos');
    case 'Delta'
        val = a.delta;
    otherwise
        val = get(mirdata(a),propName);
end