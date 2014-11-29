function demo8classification
% To get familiar with different approaches of classification using 
% MIRtoolbox, and to assess their performances. 

% Part 1. The aim of this experiment is to categorize a set of very short 
% musical excerpts according to their genres, through a supervised learning.

try
    cd classif
catch
    error('Please change current directory to ''MIRtoolboxDemos'' directory')
end

mirverbose(0)

% 1.4. Compute the mel-frequency cepstrum coefficient for each different 
% audio file
mfcc = mirmfcc('Folders');

% 1.5.  Estimate the label (i.e., genre) of each file based on a prior 
% learning using the other files of the set. Use for this purpose
% the mirclassify function.
help mirclassify

% Let's first try a classification based of mfcc, for instance, using the
% minimum distance strategy:
mirclassify(mfcc)

% The results indicates the outcomes and the total correct classification
% rate (CCR).

% 1.6. Let's try a k-nearest-neighbour strategy. For instance, for k = 5:
mirclassify(mfcc,5)

% 1.7. Use a Gaussian mixture modelling with one gaussian per class:
mirclassify(mfcc,'GMM',1)

% try also with three Gaussians per class.
mirclassify(mfcc,'GMM',3)

% As this strategy is stochastic, the results vary for every trial.
mirclassify(mfcc,'GMM',1)
mirclassify(mfcc,'GMM',1)
mirclassify(mfcc,'GMM',3)
mirclassify(mfcc,'GMM',3)

% 1.8. Carry out the classification using other features such as spectral
% centroid:
centroid = mircentroid('Folders');
mirclassify(centroid,'GMM',1)
mirclassify(centroid,'GMM',1)
mirclassify(centroid,'GMM',3)
mirclassify(centroid,'GMM',3)

% try also spectral entropy and spectral irregularity. 
entropy = mirentropy('Folders');
mirclassify(entropy,'GMM',1)
mirclassify(entropy,'GMM',1)
mirclassify(entropy,'GMM',3)
mirclassify(entropy,'GMM',3)

irregularity = mirregularity('Folders','Contrast',.1);
mirclassify(irregularity,'GMM',1)
mirclassify(irregularity,'GMM',1)
mirclassify(irregularity,'GMM',3)
mirclassify(irregularity,'GMM',3)

% Try classification based on a set of features such as:
mirclassify({entropy,centroid},'GMM',1)
mirclassify({entropy,centroid},'GMM',1)
mirclassify({entropy,centroid},'GMM',3)
mirclassify({entropy,centroid},'GMM',3)

% 1.9. By varying the features used for classification, the strategies and
% their parameters, try to find an optimal strategy that give best correct
% classification rate.
bright = mirbrightness('Folders');
rolloff = mirrolloff('Folders');
spread = mirspread('Folders');
mirclassify({bright,rolloff,spread},'GMM',3)
skew = mirskewness('Folders');
kurtosis = mirkurtosis('Folders');
flat = mirflatness('Folders');
mirclassify({skew,kurtosis,flat},'GMM',3)
for i = 1:3
     mirclassify({mfcc,centroid,skew,kurtosis,flat,entropy,irregularity,...
                 bright,rolloff,spread},'GMM',3)
end

cd ..

%%
% Part 2. In this second experiment, we will try to cluster the segments of
% an audio file according to their mutual similarity.

% 2.1.  To simplify the computation, downsample
% the audio file to 11025 Hz.
a = miraudio('czardas.wav','Sampling',11025);

% 2.2. Decompose the file into successive frames of 2 seconds with half-
% overlapping.
f = mirframe(a,2,.1);

% 2.3. Segment the file based on the novelty of the key strengths.
n = mirnovelty(mirkeystrength(f),'KernelSize',5)
p = mirpeaks(n)
s = mirsegment(a,p)

% 2.4. Compute the key strengths of each segment.
ks = mirkeystrength(s)

% 2.5. Cluster the segments according to their key strengths.
help mircluster
mircluster(s,ks)

% The k means algorithm used in the clustering is stochastic, and its
% results may vary at each run. By default, the algorithm is run 5 times 
% and the best result is selected. Try the analysis with a higher number of
% runs:
mircluster(s,ks,'Runs',10)
