function demo2timbre

try
    cd('train_set')
catch
end

a = miraudio('Folder');

z = mirzerocross(a)
disp('Ascending order of zero crossing...')
mirplay(a,'Increasing',z,'Every',5)

l = mirlowenergy(a)
disp('Ascending order of low energy rate...')
mirplay(a,'Increasing',l,'Every',5)

c = mircentroid(a)
disp('Ascending order of spectral centroid...')
mirplay(a,'Increasing',c,'Every',5)

r = mirrolloff(a)
disp('Ascending order of spectral roll-off...')
mirplay(a,'Increasing',r,'Every',5)

i = mirregularity(a)
disp('Ascending order of spectral irregularity...')
mirplay(a,'Increasing',i,'Every',5)

e = mirentropy(a)
disp('Ascending order of spectral entropy...')
mirplay(a,'Increasing',e,'Every',5)