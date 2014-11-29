help mirtoolbox
help miraudio
a = miraudio('ragtime.wav','Center','Sampling',11025)
mirplay(a)
mirspectrum(a)
mirgetdata(a)
get(a,'Sampling')

a = miraudio('laksin.wav')
a = miraudio(a,'Trim')
a = miraudio(a,'Extract',0,1)
e = mirenvelope(a)
e = mirenvelope(a,'Diff')
e = mirenvelope(a,'HalfwaveDiff')

pause, close all

s = mirspectrum('Amin3.wav')
s = mirspectrum(s,'Max',3000)
s = mirspectrum(s,'dB')
s = mirspectrum('Amin3.wav','Mel')
mirpeaks(s)

pause, close all

ac = mirautocor('Amin3.wav')
ac = mirautocor(ac,'Enhanced',2:10)
ac = mirautocor(ac,'Freq')

pause, close all

f = mirframe('ragtime.wav',.1,.5)
s = mirspectrum(f)
s = mirspectrum('ragtime.wav','frame',.1,.5)
mirpeaks(s)
mirflux(s)
mirflux(mirautocor(a,'Frame'))

pause, close all

fb = mirfilterbank('ragtime.wav','NbChannels',5)
e = mirenvelope(fb)
ae = mirautocor(e)
sa = mirsum(ae)

pause, close all

r1 = mirrms('movie1.wav','Frame')
r2 = mirrms('movie2.wav','Frame')
mirlowenergy(r1)
mirlowenergy(r2)

pause, close all

[t,ac] = mirtempo('ragtime.wav')
[t,ac] = mirtempo('czardas.wav','frame')

pause, close all

at = mirattacks('ragtime.wav')
mirattackslope(at)
mirbrightness('ragtime.wav','frame')
mirroughness('ragtime.wav','frame')

pause, close all

[p,a] = mirpitch('ragtime.wav','frame')
mirchromagram('ragtime.wav','wrap','no')
mirchromagram('ragtime.wav')
mirkey('ragtime.wav')
k = mirkey('ragtime.wav','frame')
mirmode('ragtime.wav')
m = mirmode('ragtime.wav','frame')

pause, close all

mirstat(m)
mirhisto(m)
mirexport('result.txt',k,m)