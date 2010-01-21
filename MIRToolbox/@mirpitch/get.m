function val = get(p, propName)
% GET Get properties from the MIRpitch object
% and return the value

switch propName
    case 'Amplitude'
        val = p.amplitude;
    otherwise
        val = get(mirscalar(p),propName);
end