function val = get(a, propName)
% GET Get properties from the MIRchromagram object
% and return the value

switch propName
    case 'Magnitude'
        val = get(mirdata(a),'Data');
    case 'Chroma'
        val = get(mirdata(a),'Pos');
    case 'ChromaClass'
        val = a.chromaclass;
    case 'ChromaFreq'
        val = a.chromafreq;
    case 'Register'
        val = a.register;
    case 'PitchLabel'
        val = a.plabel;
    case 'Wrap'
        val = a.wrap;
    otherwise
        val = get(mirdata(a),propName);
end