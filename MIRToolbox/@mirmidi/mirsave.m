function mirsave(m)

nmat = get(m,'Data');
n = get(m,'Name');
for k = 1:length(nmat)
    writemidi(nmat{k},[n{k},'.mid'],120,60);
end