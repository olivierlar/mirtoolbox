help mirtoolbox
help miraudio


a = miraudio('ragtime','Center','Sampling',11025,'Normal')
mirplay(a)
a = miraudio('ragtime','Extract',0,1)
mirplay(a)
miraudio('ragtime','Trim')
a1 = miraudio('pianoA4');
a2 = miraudio('pianoF4');
a3 = a1+a2;
mirplay(a3)
mirsave(a3)

f = mirframe('ragtime',1,.5)
mirplay(f)

mirenvelope('ragtime')
mirenvelope('ragtime','Tau',.05)
mirenvelope('ragtime','Diff')
mirenvelope('ragtime','HalfwaveDiff')

s = mirspectrum('pianoF4')
mirspectrum(s,'Max',3000)
mirspectrum('pianoF4','dB')
mirspectrum('pianoF4','Mel')
mirspectrum('trumpet')
mirspectrum('trumpet','Prod',2:6)

c = mircepstrum('pianoA4')
mircepstrum(c,'Freq')

mirautocor('trumpet')
ac = mirautocor('Amin3','Freq')
mirautocor(ac,'Halfwave')
mirautocor(ac,'Enhanced')
mirautocor(ac,'Enhanced',2:10)

as = mirautocor(mirspectrum('Amin3'))
ac = mirautocor('Amin3','Freq')
cp = mircepstrum('Amin3','Freq')
ac*as
ac*cp
as*cp

mirspectrum('ragtime','frame')
mirflux(ans)
mircepstrum('ragtime','frame')
mirflux(ans)

fb = mirfilterbank('ragtime','Gammatone')
mirsum(fb)
s = mirspectrum(fb)
mirsummary(s)
mirauditory('ragtime')
mirauditory('ragtime','Filterbank',20)

mirpeaks(mirspectrum('ragtime','mel'))
mirpeaks(mirspectrum('ragtime','mel','frame'),'total',1)


r1 = mirrms('movie1','Frame')
r2 = mirrms('movie2','Frame')
mirlowenergy(r1)
mirlowenergy(r2)

s = mirspectrum('ragtime','Frame',.023,.5,'Mel', 'dB')
s2 = mirspectrum(s,'AlongBands','Max',10,'Window', 0,'Resonance', 'Fluctuation')
mirsum(s2)

mironsets('ragtime')
mironsets('ragtime','Detect',0)
mironsets('ragtime','Diffenvelope')
mironsets('ragtime','diffenvelope','Contrast',.1)
mironsets('ragtime','SpectralFlux')
mironsets('ragtime','SpectralFlux','Inc','off')
mironsets('ragtime','SpectralFlux','Complex')

[t,a] = mirtempo('ragtime')
[t,a] = mirtempo('ragtime','spectrum')
[t,a] = mirtempo('ragtime','frame')

%[p s] = mirpulseclarity('ragtime')

mirattacks('ragtime')
mirattacktime('ragtime')
mirattackslope('ragtime')

t = mirtempo('czardas','frame')
st = mirstat(t)
h = mirhisto(t)
mirexport('result.txt',t)