function val = get(m, propName)
% GET Get properties from the MIRmidi object
% and return the value

switch propName
    case 'Data'
        val = m.data;
    otherwise
        val = get(mirdata(m),propName);
end