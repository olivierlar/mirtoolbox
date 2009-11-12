function v = get(a,varargin)
% GET Get properties from the MIRstruct object and return the value

switch varargin{1}
    case 'Fields'
        v = a.fields;
    case 'Data'
        v = a.data;
    case 'Tmp'
        v = a.tmp;
    case 'Stat'
        v = a.stat;
    otherwise
        v = get(mirdesign(a),varargin{:});
end