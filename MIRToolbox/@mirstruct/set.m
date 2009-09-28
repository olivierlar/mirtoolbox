function ss = set(s,varargin)
% SET Set properties from the MIRstruct object and return the value

propertyArgIn = varargin;
f = s.fields;
d = s.data;
t = s.tmp;
des = mirdesign(s);
while length(propertyArgIn) >= 2,
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch prop
        case 'Fields'
            f = val;
        case 'Data'
            d = val;
        case 'Tmp'
            t = val;
        otherwise
            des = set(des,prop,val);
    end
end
ss.fields = f;
ss.data = d;
ss.tmp = t;
ss = class(ss,'mirstruct',des);