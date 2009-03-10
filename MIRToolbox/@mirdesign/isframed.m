function b = isframed(d)

if isstruct(d.frame)
    b = 1;
else
    b = 0;
end