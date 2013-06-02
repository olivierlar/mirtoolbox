function val = get(e, propName)
% GET Get properties from the MIRenvelope object
% and return the value

switch propName
    case 'DownSampling'
        val = e.downsampl;
    case 'Halfwave'
        val = e.hwr;
    case 'Diff'
        val = e.diff;
    case 'Log';
        val = e.log;
    case 'Phase'
        val = e.phase;
    otherwise
        val = get(mirtemporal(e),propName);
end