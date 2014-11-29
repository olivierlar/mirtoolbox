function demo1pitch

% Demo 1:Pitch extraction		
% Aims: To study the diverse tools available in the toolbox for pitch
%extraction.
% In the following exercises, we will try to extract the pitch contents of various audio files using diverse techniques. All the audio files are included in the ?MIRtoolbox_Demos? folder.

%% 1. Analyis of a single chord

% Load the audio file 'Amin3.wav'.
a = miraudio('Amin3.wav')
%and play it:
mirplay(a)

% Observe the periodicities contained in the signal by computing the autocorrelation function:
ac = mirautocor(a)
% The autocorrelation function can also be displayed in the frequency domain (?Freq?), and the frequency range can be specified, for instance between 75 and 2400 Hz:
ac = mirautocor(ac, 'Freq','Min',75,'Hz','Max',2400,'Hz')

% The peaks of the curve indicates the most important periodicities:
pac = mirpeaks(ac)
% But the two peaks at the start and end of the curves are meaningless, so should be removed:
pac = mirpeaks(ac, 'NoBegin','NoEnd')

% The corresponding pitch height are given by:
p = mirpitch(pac)

% The actual numerical data (here, the frequency of each pitch) is given by:
mirgetdata(p)

% The resulting pitches can be finally played:
mirplay(p)

% So far, the results do not seem very convincing, do they?

pause, close all

%% 2. Improving the analysis

% Clearly, we need to refine the methods in order to retrieve the notes played in the chord. First of all, the autocorrelation function can be ?generalized?, i.e., ?compressed? by a factor of, for instance, 0.5 (cf. MIRtoolbox User?s Manual at the mirautocor section for an explanation of the 'Compres' option):
ac = mirautocor(a,'Compres',.5)

% The rest of the computation can actually be performed using the shortcut:
[p pac] = mirpitch(ac)  

% Look and listen to the results. 
mirgetdata(p)
mirplay(p)

% Does it sound better? What is the problem this time?

pause, close all

%% 3. Improving the analysis, again

% In fact, the autocorrelation function contains a lot of harmonics that need to be removed. For that purpose, the autocorrelation function can be ?enhanced? (cf. MIRtoolbox User?s Manual at the mirautocor section for an explanation of the 'Enhanced' option):
ac = mirautocor(ac,'Enhanced')

% Carry out the rest of the computation as in section 2.2. 
[p pac] = mirpitch(ac)
mirgetdata(p)
mirplay(p)

% What do you think of the results?

pause, close all

%% 4. Improving the analysis, still

% An additional improvement consists in first decomposing the audio into two channels using the command:
fb = mirfilterbank(a,'2Channels') 

% Compute the autocorrelation on each channel:
ac = mirautocor(fb,'Compres',.5)

% The autocorrelation of each channel can then be summed together:
ac = mirsum(ac)

% And the enhancement can be performed afterwards:
ac = mirautocor(ac,'Enhanced')

% And the rest of the computation follows the same principle than before.
[p ac] = mirpitch(ac)

% The result should be better this way.
mirgetdata(p)
mirplay(p)

% Hopefully, the whole chain of operation can be performed automatically using the simple command:
[p ac] = mirpitch(a)
[p ac] = mirpitch(a,'Cepstrum')
[p ac] = mirpitch(a,'AutocorSpectrum')

pause, close all

%% 5. Frame-based analysis

% Let?s analyze a longer musical sequence, such as ragtime.wav:
a = miraudio('ragtime.wav');
mirplay(a)

% A direct analysis of the recording using the previous command
[p ac] = mirpitch(a) 
% is not concluding, as it tries to find the pitches contained in the whole signal.
mirgetdata(p)
mirplay(p)

% It is better to estimate the pitch content frame-by-frame:
[p ac] = mirpitch(a,'Frame')
mirgetdata(p)

% Don?t forget to listen to the results.
mirplay(p)

pause, close all

%% 6. Segment-based analysis

% But it might get more sense to first decompose the signal into notes, using the commands:
o = mironsets(a,'Attacks')
% (The 'Attacks' option enables to segment at the beginning of each note
% attack phase, i.e., where the energy is minimum, instead of segmenting where the energy is maximum.)
mirplay(o)

sg = mirsegment(a,o)
mirplay(sg)

% and to compute the pitch on each successive segment:
[p ac] = mirpitch(sg)
mirgetdata(p)

% Again, don?t forget to listen to the results.
mirplay(p)

pause, close all

%% 7. Monody analysis

% Conclude with the analysis of a Finnish folk song, Läksin minä kesäyönä, whose excerpt is recorded in the file ?laksin.wav?. As the melody is sung by one voice only, you can use the ?Mono? option in mirpitch.
mirplay('laksin.wav')
o = mironsets('laksin.wav','Attacks','Contrast',.1)
mirplay(o)
sg = mirsegment('laksin.wav',o)
mirplay(sg)
p = mirpitch(sg,'mono')
mirgetdata(p)
mirplay(p)