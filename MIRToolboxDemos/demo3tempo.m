function demo3tempo
% To get familiar with tempo estimation from audio using the MIR Toolbox.
% To assess the performance of the tempo estimation method.

% 1. Let's investigate the different stages needed for tempo estimation. 
 
d = miraudio('ragtime.wav')
mirenvelope(d)
e = mirenvelope(d,'Halfwavediff')

% Decompose the audio file with a filter bank.
f = mirfilterbank(d)

% Calculate also a half-wave rectified differentiated envelope.
ee = mirenvelope(f,'HalfwaveDiff')

% Sum the frequency channels.ok
s = mirsum(ee,'Centered') 

d2 = miraudio('vivaldi.wav')
f2 = mirfilterbank(d2)
ee2 = mirenvelope(f2,'HalfwaveDiff')
s2 = mirsum(ee2,'Centered') 

% Calculate the autocorrelation function.
ac = mirautocor(s) 
 
% Apply the resonance model to the autocorrelation function.
ac = mirautocor(s,'Resonance') 
 
% Find peaks in the autocorrelation function.
p = mirpeaks(ac) 
mirgetdata(p)

% Get the period of the peaks.
t = mirtempo(p,'Total',1)

display('Strike any key to continue...');
pause
close all

% 2. All the functions we used are integrated into the function tempo.
help mirtempo

% For instance, we can simply write:
[t,ac] = mirtempo('ragtime.wav')

% As you can see in the help, the resonance is integrated by default in the
% tempo function. To toggle off the use of the resonance function, type:
[t,ac] = mirtempo('ragtime.wav','Resonance',0)

[t,ac] = mirtempo('ragtime.wav','total',5)

display('Strike any key to continue...');
pause
close all

% 3. The excerpt 'laksin' and 'czardas' have variable tempi. Use frame-
% based tempo analysis to estimate the variation of the tempi. 
% Apply to this end the tempo command with the 'frame' option. 
[t1,p1] = mirtempo('laksin.wav','Frame')
[t2,p2] = mirtempo('czardas.wav','Frame')

% What is the range of variation of tempi?
help mirhisto
h1 = mirhisto(t1)
h2 = mirhisto(t2)