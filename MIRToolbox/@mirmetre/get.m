function val = get(m, propName)
% GET Get properties from the MIRmetre object
% and return the value

switch propName
    case 'Meters'
        val = m.meters;
end