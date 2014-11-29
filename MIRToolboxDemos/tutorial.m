help mirtoolbox
help miraudio


a = miraudio('ragtime.wav','Center','Sampling',11025,'Normal')
mirplay(a)
a = miraudio('ragtime.wav','Extract',0,1)
mirplay(a)
miraudio('ragtime.wav','Trim')
a1 = miraudio('pianoA4.wav');
a2 = miraudio('pianoF4.wav');
a3 = a1+a2;
mirplay(a3)
mirsave(a3)

f = mirframe('ragtime.wav',1,.5)
mirplay(f)

mirenvelope('ragtime.wav')
mirenvelope('ragtime.wav','Tau',.05)
mirenvelope('ragtime.wav','Diff')
mirenvelope('ragtime.wav','HalfwaveDiff')

s = mirspectrum('pianoF4.wav')
mirspectrum(s,'Max',3000)
mirspectrum('pianoF4.wav','dB')
mirspectrum('pianoF4.wav','Mel')
mirspectrum('trumpet.wav')
mirspectrum('trumpet.wav','Prod',2:6)

c = mircepstrum('pianoA4.wav')
mircepstrum(c,'Freq')

mirautocor('trumpet.wav')
ac = mirautocor('Amin3.wav','Freq')
mirautocor(ac,'Halfwave')
mirautocor(ac,'Enhanced')
mirautocor(ac,'Enhanced',2:10)

as = mirautocor(mirspectrum('Amin3.wav'))
ac = mirautocor('Amin3.wav','Freq')
cp = mircepstrum('Amin3.wav','Freq')
ac*as
ac*cp
as*cp

mirspectrum('ragtime.wav','frame')
mirflux(ans)
mircepstrum('ragtime.wav','frame')
mirflux(ans)

fb = mirfilterbank('ragtime.wav','Gammatone')
mirsum(fb)
s = mirspectrum(fb)
mirsummary(s)
mirauditory('ragtime.wav')
mirauditory('ragtime.wav','Filterbank',20)

close all

mirpeaks(mirspectrum('ragtime.wav','mel'))
mirpeaks(mirspectrum('ragtime.wav','mel','frame'),'total',1)


r1 = mirrms('movie1.wav','Frame')
r2 = mirrms('movie2.wav','Frame')
mirlowenergy(r1)
mirlowenergy(r2)

s = mirspectrum('ragtime.wav','Frame',.023,.5,'Mel', 'dB')
s2 = mirspectrum(s,'AlongBands','Max',10,'Window', 0,'Resonance', 'Fluctuation')
mirsum(s2)

mironsets('ragtime.wav')
mironsets('ragtime.wav','Detect',0)
mironsets('ragtime.wav','Diffenvelope')
mironsets('ragtime.wav','diffenvelope','Contrast',.1)
mironsets('ragtime.wav','SpectralFlux')
mironsets('ragtime.wav','SpectralFlux','Inc','off')
mironsets('ragtime.wav','SpectralFlux','Complex')

[t,a] = mirtempo('ragtime.wav')
[t,a] = mirtempo('ragtime.wav','spectrum')
[t,a] = mirtempo('ragtime.wav','frame')

%[p s] = mirpulseclarity('ragtime')

mirattacks('ragtime.wav')
mirattacktime('ragtime.wav')
mirattackslope('ragtime.wav')

t = mirtempo('czardas.wav','frame')
st = mirstat(t)
h = mirhisto(t)
mirexport('result.txt',t)