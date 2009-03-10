function val = get(a, propName)
% GET Get properties from the MIRclassify object
% and return the value

switch propName
    case 'Correct'
        val = a.correct;
    case 'Data'
        val = a.correct;
    case 'Classes'
        val = a.classes;
    otherwise
        error([propName,' is not a valid MIRdata property'])
end