help mirtoolbox
help miraudio
a = miraudio('ragtime','Center','Sampling',11025)
mirplay(a)
mirspectrum(a)
mirgetdata(a)
get(a,'Sampling')

a = miraudio('laksin')
a = miraudio(a,'Trim')
a = miraudio(a,'Extract',0,1)
e = mirenvelope(a)
e = mirenvelope(a,'Diff')
e = mirenvelope(a,'HalfwaveDiff')

pause, close all

s = mirspectrum('Amin3')
s = mirspectrum(s,'Max',3000)
s = mirspectrum(s,'dB')
s = mirspectrum('Amin3','Mel')
mirpeaks(s)

pause, close all

ac = mirautocor('Amin3')
ac = mirautocor(ac,'Enhanced',2:10)
ac = mirautocor(ac,'Freq')

pause, close all

f = mirframe('ragtime',.1,.5)
s = mirspectrum(f)
s = mirspectrum('ragtime','frame',.1,.5)
mirpeaks(s)
mirflux(s)
mirflux(mirautocor(a,'Frame'))

pause, close all

fb = mirfilterbank('ragtime','NbChannels',5)
e = mirenvelope(fb)
ae = mirautocor(e)
sa = mirsum(ae)

pause, close all

r1 = mirrms('movie1','Frame')
r2 = mirrms('movie2','Frame')
mirlowenergy(r1)
mirlowenergy(r2)

pause, close all

[t,ac] = mirtempo('ragtime')
[t,ac] = mirtempo('czardas','frame')

pause, close all

at = mirattacks('ragtime')
mirattackslope(at)
mirbrightness('ragtime','frame')
mirroughness('ragtime','frame')

pause, close all

[p,a] = mirpitch('ragtime','frame')
mirchromagram('ragtime','wrap','no')
mirchromagram('ragtime')
mirkey('ragtime')
k = mirkey('ragtime','frame')
mirmode('ragtime')
m = mirmode('ragtime','frame')

pause, close all

mirstat(m)
mirhisto(m)
mirexport('result.txt',k,m)