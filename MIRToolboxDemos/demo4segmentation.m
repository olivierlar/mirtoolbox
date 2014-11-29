function demo4segmentation
% To get familiar with some approaches of segmentation of audio files 
% using MIRtoolbox. 
 
% 1. Load an audio file (for instance, guitar.wav). 
a = miraudio('guitar.wav');

% 2. We will perform the segmentation strategy as proposed in (Foote &
% Cooper, 2003). First, decompose the file into successive frames of 50 ms 
% without overlap. 
help mirframe
fr = mirframe(a,0.05,1)

% 3. Compute the spectrum representation (FFT) of the frames.
sp = mirspectrum(fr)
clear fr
% (Remove from the memory any data that will not be used any more.)

% 4. Compute the similarity matrix that shows the similarity between the 
% spectrum of different frames.
help mirsimatrix
sm = mirsimatrix(sp)
clear sp
% Look at the structures shown in the matrix and find the relation with the
% structure heard when listening to the extract. 

% 5. Estimate the novelty score related to the similarity matrix. It 
% consists in a convolution of the diagonal of the matrix with a
% checker-board Gaussian kernel. Use the novelty function for that purpose. 
help mirnovelty
nv = mirnovelty(sm)

% 6. Detect the peaks in the novelty score.
help mirpeaks
p1 = mirpeaks(nv)

% You can change the threshold value of the peak picker function in order to 
% get better results. 
p2 = mirpeaks(nv,'Contrast',0.01)

clear nv

% 7. Segment the original audio file using the peaks as position for
% segmentation. 
help mirsegment
s1 = mirsegment(a,p1)
clear p1

% 8. Listen to the results.
mirplay(s1)

%s2 = mirsegment(a,p2)
%clear p2
%mirplay(s2)

% 9. Compute the similarity matrix of this obtained segmentation, in order
% to view the relationships between the different segments and their
% possible clustering into higher-level groups. 
mirsimatrix(s1,'Similarity')
clear s1
%mirsimatrix(s2)
%clear s2

display('Strike any key to continue...');
pause
close all

% 10. Change the size of the kernel used in the novelty function, in order
% to obtain segmentations of different levels of detail, from detailed
% analysis of the local texture, to very simple segmentation of the whole 
% piece.
n100 = mirnovelty(sm,'KernelSize',100)
n50 = mirnovelty(sm,'KernelSize',50)
n10 = mirnovelty(sm,'KernelSize',10)
clear sm
% As you can see, the smaller the gaussian kernel is, the more peaks can be
% found in the novelty score. Indeed, if the kernel is small, the cumulative
% multiplication of its elements with the superposed elements in the
% similarity matrix may vary more easily, throughout the progressive
% sliding of the kernel along the diagonal of the similarity matrix, and
% local change of texture may be more easily detected. On the contrary,
% when the kernel is large, only large-scale change of texture are
% detected.

display('Strike any key to continue...');
pause
close all

p100 = mirpeaks(n100,'NoBegin','NoEnd')
clear n100
p50 = mirpeaks(n50,'NoBegin','NoEnd')
clear n50
p10 = mirpeaks(n10,'NoBegin','NoEnd')
clear n10
s100 = mirsegment(a,p100)
clear p100
mirplay(s100)
clear s100
s50 = mirsegment(a,p50)
clear p50
mirplay(s50)
clear s50
s10 = mirsegment(a,p10)
clear p10
mirplay(s10)
clear s10

display('Strike any key to continue...');
pause
close all

% One more compact way of writing these commands is as follows:
mirsegment(a,'Novelty')
mirsegment(a,'Novelty','Contrast',0.01)
mirsegment(a,'Novelty','KernelSize',100)

display('Strike any key to continue...');
pause
close all

% Besides, if you want to see the novelty curve with the peaks, just add a
% second output:
[s50 p50] = mirsegment(a,'Novelty','KernelSize',50)
clear s50 p50
[s10 p10] = mirsegment(a,'Novelty','KernelSize',10)
clear a s10 p10

display('Strike any key to continue...');
pause
close all


% 11. Try the whole process with MFCC instead of spectrum analysis. Take the
% first ten MFCC for instance.
help mirsegment
% The segment function can simply be called as follows:
sc = mirsegment('czardas.wav','Novelty','MFCC','Rank',1:10)
clear sc

% Here are some other examples of use:
[ssp p m b] = mirsegment('valse_triste_happy.wav','Spectrum',...
                                'KernelSize',150,'Contrast',.1)
clear p m b
mirplay(ssp)
clear ssp

display('Strike any key to continue...');
pause
close all

[smfcc2 p m a] = mirsegment('valse_triste_happy.wav','MFCC',2:10,...
                                'KernelSize',150,'Contrast',.1)
clear p m a
mirplay(smfcc2)