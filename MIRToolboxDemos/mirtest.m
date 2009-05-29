function mirtest(audio)

if not(nargin)
    audio = 'ragtime';
end

%% Version 1.0

mirpeaks(mirspectrum(audio,'Mel'))
f = mirframe(audio,.5,.5)
mirpeaks(mirspectrum(f,'Mel'),'Total',1)
mirpeaks(mirautocor(f),'Total',1,'NoBegin')
mirpeaks(mirspectrum(f),'Total',1,'NoBegin')
mirpeaks(mirchromagram(f),'Total',1)
mirpeaks(mirkeystrength(f),'Total',1)
%mirpeaks(mirfluctuation(f),'Total',1,'NoBegin') %Not implemented yet..
[a,b,c] = mirkey(audio)

%%
pause
clear a b c f
close all

f = mirfeatures(audio);
sf = mirstat(f);
f.dynamics.rms{1}
f.rhythm.fluctuation.peak{1}
f.rhythm.fluctuation.centroid{1}
f.rhythm.tempo{:}
f.rhythm.attack.time{:}
f.rhythm.attack.slope{1}
sf.dynamics.rms
sf.rhythm.fluctuation.peak
sf.rhythm.fluctuation.centroid
sf.rhythm.tempo
sf.rhythm.attack.time
sf.rhythm.attack.slope

pause
close all

f.timbre.zerocross{:}
f.timbre.centroid{:}
f.timbre.brightness{:}
f.timbre.spread{:}
f.timbre.skewness{:}
f.timbre.kurtosis{:}
f.timbre.rolloff95{:}
f.timbre.rolloff85{:}
f.timbre.spectentropy{:}
f.timbre.flatness{:}
sf.timbre.zerocross
sf.timbre.centroid
sf.timbre.brightness
sf.timbre.spread
sf.timbre.skewness
sf.timbre.kurtosis
sf.timbre.rolloff95
sf.timbre.rolloff85
sf.timbre.spectentropy
sf.timbre.flatness

pause
close all

f.timbre.roughness{:}
f.timbre.irregularity{:}
f.timbre.inharmonicity{:}
f.timbre.mfcc{:}
f.timbre.dmfcc{:}
f.timbre.ddmfcc{:}
f.timbre.lowenergy{:}
sf.timbre.roughness
sf.timbre.irregularity
sf.timbre.inharmonicity
sf.timbre.mfcc
sf.timbre.dmfcc
sf.timbre.ddmfcc
sf.timbre.lowenergy

pause
close all

f.timbre.spectralflux{:}
f.pitch.salient{:}
f.pitch.chromagram.peak{:}
f.pitch.chromagram.centroid{:}
f.tonal.keyclarity{:}
f.tonal.mode{:}
f.tonal.hcdf{:}
sf.timbre.spectralflux
sf.pitch.salient
sf.pitch.chromagram.peak
sf.pitch.chromagram.centroid
sf.tonal.keyclarity
sf.tonal.mode
sf.tonal.hcdf

mirexport('resultdemo.txt',sf)
mirexport('resultdemo.arff',f)

%% Version 1.1

pause
clear f
close all

mirlength(audio)
s = mirspectrum(audio,'cents','Min',50)
s = mirspectrum(s,'Collapsed')
mirspectrum(s,'Gauss')
ss = mirspectrum(s,'Smooth')
p = mirpeaks(ss,'Extract')
mirkurtosis(p)
[le,f] = mirlowenergy(audio,'ASR')
p = mirpitch(audio,'frame')
mirpitch(p,'median')
mirauditory(audio)
mirroughness('ragtime')

%%
pause
clear s ss p le f
close all

fb = mirfilterbank('Design','NbChannels',5)
f = mirfeatures(fb);
%sf = mirstat(f);
f = mireval(f,audio)

f.dynamics.rms
f.rhythm.fluctuation.peak
f.rhythm.fluctuation.centroid
f.rhythm.tempo
f.rhythm.attack.time
f.rhythm.attack.slope
%sf.dynamics.rms
%sf.rhythm.fluctuation.peak
%sf.rhythm.fluctuation.centroid
%sf.rhythm.tempo
%sf.rhythm.attack.time
%sf.rhythm.attack.slope

pause
close all

f.timbre.zerocross
f.timbre.centroid
f.timbre.brightness
f.timbre.spread
f.timbre.skewness
f.timbre.kurtosis
f.timbre.rolloff95
f.timbre.rolloff85
f.timbre.spectentropy
%sf.timbre.zerocross
%sf.timbre.centroid
%sf.timbre.brightness
%sf.timbre.spread
%sf.timbre.skewness
%sf.timbre.kurtosis
%sf.timbre.rolloff95
%sf.timbre.rolloff85
%sf.timbre.spectentropy

pause
close all

f.timbre.flatness
f.timbre.roughness
f.timbre.irregularity
f.timbre.inharmonicity
f.timbre.mfcc
f.timbre.dmfcc
f.timbre.ddmfcc
f.timbre.lowenergy
%sf.timbre.flatness
%sf.timbre.roughness
%sf.timbre.irregularity
%sf.timbre.inharmonicity
%sf.timbre.mfcc
%sf.timbre.dmfcc
%sf.timbre.ddmfcc
%sf.timbre.lowenergy

pause
close all

f.timbre.spectralflux
f.pitch.salient
f.pitch.chromagram.peak
f.pitch.chromagram.centroid
f.tonal.keyclarity
f.tonal.mode
f.tonal.hcdf
%sf.timbre.spectralflux
%sf.pitch.salient
%sf.pitch.chromagram.peak
%sf.pitch.chromagram.centroid
%sf.tonal.keyclarity
%sf.tonal.mode
%sf.tonal.hcdf

pause
close all

f.tmp.fluctuation
f.tmp.onsets
f.tmp.attacks
f.tmp.s
f.tmp.pitch
f.tmp.mfcc
f.tmp.dmfcc
f.tmp.chromagram
f.tmp.keystrengths
%sf.tmp.fluctuation
%sf.tmp.onsets
%sf.tmp.attacks
%sf.tmp.s
%sf.tmp.pitch
%sf.tmp.mfcc
%sf.tmp.dmfcc
%sf.tmp.chromagram
%sf.tmp.keystrengths

%%

pause
clear s f
close all

[s,n] = mirsegment(audio,'KernelSize',64)
f = mirfeatures(s)
%sf = mirstat(f);

f.dynamics.rms
f.rhythm.fluctuation.peak
f.rhythm.fluctuation.centroid
f.rhythm.tempo
f.rhythm.attack.time
f.rhythm.attack.slope
%sf.dynamics.rms
%sf.rhythm.fluctuation.peak
%sf.rhythm.fluctuation.centroid
%sf.rhythm.tempo
%sf.rhythm.attack.time
%sf.rhythm.attack.slope

pause
close all

f.timbre.zerocross
f.timbre.centroid
f.timbre.brightness
f.timbre.spread
f.timbre.skewness
f.timbre.kurtosis
f.timbre.rolloff95
f.timbre.rolloff85
f.timbre.spectentropy
%sf.timbre.zerocross
%sf.timbre.centroid
%sf.timbre.brightness
%sf.timbre.spread
%sf.timbre.skewness
%sf.timbre.kurtosis
%sf.timbre.rolloff95
%sf.timbre.rolloff85
%sf.timbre.spectentropy

pause
close all

f.timbre.flatness
f.timbre.roughness
f.timbre.irregularity
f.timbre.inharmonicity
f.timbre.mfcc
f.timbre.dmfcc
f.timbre.ddmfcc
f.timbre.lowenergy
%sf.timbre.flatness
%sf.timbre.roughness
%sf.timbre.irregularity
%sf.timbre.inharmonicity
%sf.timbre.mfcc
%sf.timbre.dmfcc
%sf.timbre.ddmfcc
%sf.timbre.lowenergy

pause
close all

f.timbre.spectralflux
f.pitch.salient
f.pitch.chromagram.peak
f.pitch.chromagram.centroid
f.tonal.keyclarity
f.tonal.mode
f.tonal.hcdf
%sf.timbre.spectralflux
%sf.pitch.salient
%sf.pitch.chromagram.peak
%sf.pitch.chromagram.centroid
%sf.tonal.keyclarity
%sf.tonal.mode
%sf.tonal.hcdf

pause
close all

f.tmp.fluctuation
f.tmp.onsets
f.tmp.attacks
f.tmp.s
f.tmp.pitch
f.tmp.mfcc
f.tmp.dmfcc
f.tmp.chromagram
f.tmp.keystrengths
%sf.tmp.fluctuation
%sf.tmp.onsets
%sf.tmp.attacks
%sf.tmp.s
%sf.tmp.pitch
%sf.tmp.mfcc
%sf.tmp.dmfcc
%sf.tmp.chromagram
%sf.tmp.keystrengths

%%
%pause
%clear f n
%close all
%
%fb = mirfilterbank(s,'NbChannels',5)
%f = mirfeatures(fb)