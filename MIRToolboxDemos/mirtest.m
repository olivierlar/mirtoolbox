function mirtest(audio)

if not(nargin)
    audio = 'ragtime.wav';
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
%pause
clear a b c f
close all

f = mirfeatures(audio);
sf = mirstat(f);

f.dynamics.rms{1}
f.fluctuation.peak{1}
f.fluctuation.centroid{1}
f.rhythm.tempo{:}
f.rhythm.attack.time{:}
f.rhythm.attack.slope{1}
sf.dynamics.rms
sf.fluctuation.peak
sf.fluctuation.centroid
sf.rhythm.tempo
sf.rhythm.attack.time
sf.rhythm.attack.slope

%pause
close all

f.timbre.zerocross{:}
f.spectral.centroid{:}
f.spectral.brightness{:}
f.spectral.spread{:}
f.spectral.skewness{:}
f.spectral.kurtosis{:}
f.spectral.rolloff95{:}
f.spectral.rolloff85{:}
f.spectral.spectentropy{:}
f.spectral.flatness{:}
sf.timbre.zerocross
sf.spectral.centroid
sf.spectral.brightness
sf.spectral.spread
sf.spectral.skewness
sf.spectral.kurtosis
sf.spectral.rolloff95
sf.spectral.rolloff85
sf.spectral.spectentropy
sf.spectral.flatness

%pause
close all

f.spectral.roughness{:}
f.spectral.irregularity{:}
%f.timbre.inharmonicity{:}
f.spectral.mfcc{:}
f.spectral.dmfcc{:}
f.spectral.ddmfcc{:}
f.timbre.lowenergy{:}
sf.spectral.roughness
sf.spectral.irregularity
%sf.timbre.inharmonicity
sf.spectral.mfcc
sf.spectral.dmfcc
sf.spectral.ddmfcc
sf.timbre.lowenergy

%pause
close all

f.timbre.spectralflux{:}
%f.pitch.salient{:}
f.tonal.chromagram.peak{:}
f.tonal.chromagram.centroid{:}
f.tonal.keyclarity{:}
f.tonal.mode{:}
f.tonal.hcdf{:}
sf.timbre.spectralflux
%sf.pitch.salient
sf.tonal.chromagram.peak
sf.tonal.chromagram.centroid
sf.tonal.keyclarity
sf.tonal.mode
sf.tonal.hcdf

mirexport('resultdemo.txt',sf)
mirexport('resultdemo.arff',f)

%% Version 1.1

%pause
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
mirroughness('ragtime.wav')

%%
%pause
clear s ss p le f
close all

fb = mirfilterbank('Design','NbChannels',5)
f = mirfeatures(fb);
%sf = mirstat(f);
f = mireval(f,audio)

f.dynamics.rms{1}
f.fluctuation.peak{1}
f.fluctuation.centroid{1}
f.rhythm.tempo{:}
f.rhythm.attack.time{:}
f.rhythm.attack.slope{1}

%pause
close all

f.timbre.zerocross{:}
f.spectral.centroid{:}
f.spectral.brightness{:}
f.spectral.spread{:}
f.spectral.skewness{:}
f.spectral.kurtosis{:}
f.spectral.rolloff95{:}
f.spectral.rolloff85{:}
f.spectral.spectentropy{:}
f.spectral.flatness{:}

%pause
close all

f.spectral.roughness{:}
f.spectral.irregularity{:}
%f.timbre.inharmonicity{:}
f.spectral.mfcc{:}
f.spectral.dmfcc{:}
f.spectral.ddmfcc{:}
f.timbre.lowenergy{:}

%pause
close all

f.timbre.spectralflux{:}
%f.pitch.salient{:}
f.tonal.chromagram.peak{:}
f.tonal.chromagram.centroid{:}
f.tonal.keyclarity{:}
f.tonal.mode{:}
f.tonal.hcdf{:}

%pause
close all

mirfeatures('Folder')