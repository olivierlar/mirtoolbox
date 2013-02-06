function val = get(p, propName)
% GET Get properties from the MIRpitch object
% and return the value

switch propName
    case 'Amplitude'
        val = p.amplitude;
    case 'Start'
        val = p.start;
    case 'End'
        val = p.end;
    case 'Mean'
        val = p.mean;
    case 'Degrees'
        val = p.degrees;
    case 'Stable'
        val = p.stable;
    otherwise
        val = get(mirscalar(p),propName);
end