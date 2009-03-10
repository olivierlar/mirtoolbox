function t = mirtype(x)

if isa(x,'mirdesign')
    t = get(x,'Type');
else
    t = class(x);
end