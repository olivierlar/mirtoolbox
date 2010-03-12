function b = isframed(d)

if isstruct(d.frame) && not(isfield(d.frame,'chunknow'))
    b = 1;
else
    b = 0;
end