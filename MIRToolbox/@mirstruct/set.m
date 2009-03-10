function a = set(a,varargin)
% SET Set properties from the MIRstruct object and return the value

propertyArgIn = varargin;
while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch prop
        case 'Fields'
            a.fields = val;
        case 'Data'
            a.data = val;
        case 'Tmp'
            a.tmp = val;
       otherwise
           error('Unknown MIRstruct property')
    end
end