function val = get(a, propName)
% GET Get properties from the MIRautocor object
% and return the value

switch propName
    case 'Coeff'
        val = get(mirdata(a),'Data');
    case 'Delay'
        val = get(mirdata(a),'Pos');
    case 'Lag'
        val = get(mirdata(a),'Pos');
    case 'FreqDomain'
        val = a.freq;
    case 'OfSpectrum'
        val = a.ofspectrum;
    case 'Window'
        val = a.window;
    case 'Resonance'
        val = a.resonance;
    case 'Input'
        val = a.input;
    case 'Phase'
        val = a.phase;
    otherwise
        val = get(mirdata(a),propName);
end