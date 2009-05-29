function demo8classification
% To get familiar with different approaches of classification using 
% MIRtoolbox, and to assess their performances. 

% Part 1. The aim of this experiment is to categorize a set of very short 
% musical excerpts according to their genres, through a supervised learning.

% 1.3. Select the training set for current directory.
try
    cd train_set
catch
    error('Please change current directory to ''MIRtoolboxDemos'' directory')
end

% Load all the files of the folder into one audio structure (called for 
% instance training), and associate for each folder a label defined by the
% first two letters of the respective file name.
train = miraudio('Folder','Label',1:2);
cd ..

% In the same way, select the testing set for current directory, and load
% all the files including their labels:
cd test_set
test = miraudio('Folder','Label',1:2);
cd ..

% 1.4. Compute the mel-frequency cepstrum coefficient for each different 
% audio file of both sets:
mfcc_train = mirmfcc(train);
mfcc_test = mirmfcc(test);

% 1.5.  Estimate the label (i.e., genre) of each file from the testing set,
% based on a prior learning using the training set. Use for this purpose
% the classify function.
help mirclassify

% Let's first try a classification based of mfcc, for instance, using the
% minimum distance strategy:
mirclassify(test,mfcc_test,train,mfcc_train)

% The results indicates the outcomes and the total correct classification
% rate (CCR).

% 1.6. Let's try a k-nearest-neighbour strategy. For instance, for k = 5:
mirclassify(test,mfcc_test,train,mfcc_train,5)

% 1.7. Use a Gaussian mixture modelling with one gaussian per class:
mirclassify(test,mfcc_test,train,mfcc_train,'GMM',1)

% try also with three Gaussians per class.
mirclassify(test,mfcc_test,train,mfcc_train,'GMM',3)

% As this strategy is stochastic, the results vary for every trial.
mirclassify(test,mfcc_test,train,mfcc_train,'GMM',1)
mirclassify(test,mfcc_test,train,mfcc_train,'GMM',1)
mirclassify(test,mfcc_test,train,mfcc_train,'GMM',3)
mirclassify(test,mfcc_test,train,mfcc_train,'GMM',3)

% 1.8. Carry out the classification using other features such as spectral
% centroid:
spectrum_train = mirspectrum(train);
spectrum_test = mirspectrum(test);
centroid_train = mircentroid(spectrum_train);
centroid_test = mircentroid(spectrum_test);
mirclassify(test,centroid_test,train,centroid_train,'GMM',1)
mirclassify(test,centroid_test,train,centroid_train,'GMM',1)
mirclassify(test,centroid_test,train,centroid_train,'GMM',3)
mirclassify(test,centroid_test,train,centroid_train,'GMM',3)

% try also spectral entropy and spectral irregularity. 
entropy_train = mirentropy(spectrum_train);
entropy_test = mirentropy(spectrum_test);
mirclassify(test,entropy_test,train,entropy_train,'GMM',1)
mirclassify(test,entropy_test,train,entropy_train,'GMM',1)
mirclassify(test,entropy_test,train,entropy_train,'GMM',3)
mirclassify(test,entropy_test,train,entropy_train,'GMM',3)

irregularity_train = mirregularity(spectrum_train,'Contrast',.1);
irregularity_test = mirregularity(spectrum_test,'Contrast',.1);
mirclassify(test,irregularity_test,train,irregularity_train,'GMM',1)
mirclassify(test,irregularity_test,train,irregularity_train,'GMM',1)
mirclassify(test,irregularity_test,train,irregularity_train,'GMM',3)
mirclassify(test,irregularity_test,train,irregularity_train,'GMM',3)

% Try classification based on a set of features such as:
mirclassify(test,{entropy_test,centroid_test},...
         train,{entropy_train,centroid_train},'GMM',1)
mirclassify(test,{entropy_test,centroid_test},...
         train,{entropy_train,centroid_train},'GMM',1)
mirclassify(test,{entropy_test,centroid_test},...
         train,{entropy_train,centroid_train},'GMM',3)
mirclassify(test,{entropy_test,centroid_test},...
         train,{entropy_train,centroid_train},'GMM',3)

% 1.9. By varying the features used for classification, the strategies and
% their parameters, try to find an optimal strategy that give best correct
% classification rate.
bright_train = mirbrightness(spectrum_train);
bright_test = mirbrightness(spectrum_test);
rolloff_train = mirbrightness(spectrum_train);
rolloff_test = mirbrightness(spectrum_test);
spread_train = mirspread(spectrum_train);
spread_test = mirspread(spectrum_test);
mirclassify(test,{bright_test,rolloff_test,spread_test},...
         train,{bright_train,rolloff_train,spread_train},'GMM',3)
skew_train = mirskewness(spectrum_train);
skew_test = mirskewness(spectrum_test);
kurtosis_train = mirkurtosis(spectrum_train);
kurtosis_test = mirkurtosis(spectrum_test);
flat_train = mirflatness(spectrum_train);
flat_test = mirflatness(spectrum_test);
mirclassify(test,{skew_test,kurtosis_test,flat_test},...
         train,{skew_train,kurtosis_train,flat_train},'GMM',3)
for i = 1:3
     mirclassify(test,{mfcc_test,centroid_test,skew_test,kurtosis_test,...
               flat_test,entropy_test,irregularity_test,...
               bright_test,rolloff_test,spread_test},...
              train,{mfcc_train,centroid_train,skew_train,kurtosis_train,...
                    flat_train,entropy_train,irregularity_train,...
                    bright_train,rolloff_train,spread_train},'GMM',3)
end

% You can also try to change the size of the training and testing sets (by
% simply interverting them for instance). 
for i = 1:3
    mirclassify(train,{mfcc_train,centroid_train,skew_train,kurtosis_train,...
                flat_train,entropy_train,irregularity_train,...
                bright_train,rolloff_train,spread_train},...
             test,{mfcc_test,centroid_test,skew_test,kurtosis_test,...
                   flat_test,entropy_test,irregularity_test,...
                   bright_test,rolloff_test,spread_test},'GMM',3)
end

%%
% Part 2. In this second experiment, we will try to cluster the segments of
% an audio file according to their mutual similarity.

% 2.1.  To simplify the computation, downsample
% the audio file to 11025 Hz.
a = miraudio('czardas','Sampling',11025);

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
