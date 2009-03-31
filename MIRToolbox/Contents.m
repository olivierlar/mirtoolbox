% MIRtoolbox
% Version 1.1.23 31-March-2009
%
% A complete documentation is available in the downloaded folder and online.
%			http://www.jyu.fi/music/coe/materials/mirtoolbox
%
% A more detailed documentation of each function is available using the
% help command. For instance, type help miraudio.
%
%AUDIO MANIPULATION
% miraudio         - Loads and return waveform
% mirframe         - Decomposes into successive frames
% mirsegment       - Segments
% mirplay          - Plays
% mirlength        - Temporal length
%
%DATA OUTPUT
% mirgetdata       - Return any data as a structure that can be used for
%                   further computation in Matlab
% mirsave          - Save audio signals into audio files
% mirexport        - Export the analytical results to a text file
% mirstat          - Returns statistics of any feature
% mirfeatures      - Compute a large range of features
%
%TIME DOMAIN ANALYSIS
% mirrms           - Root mean square energy
% mirlowenergy     - Number of frames with lower than average energy
% mirautocor       - Autocorrelation
%
%SPECTRUM ANALYSIS
% mirspectrum      - FFT spectrum in frequency or mel-bank domains
% mirbrightness    - Spectral brightness (high-frequency rate)
% mirrolloff       - Spectral rolloff (frequency above which is located a 
%                       certain amount of energy)
%
%AUDITORY MODELLING
% mirfilterbank    - Decomposes audio signals through a bank of filters
% mirenvelope      - Signal envelope (global shape of the waveform)
% mirsum           - Sums the envelopes of a filterbank
% mirsummary       - Sums autocorrelations, spectrums, etc. of filterbanks
%
%PITCH
% mirpitch         - Pitch frequencies
% mircepstrum      - Cepstrum representation (showing periodicities)
%
%TONALITY
% mirchromagram    - Chromagram (distribution of energy along pitches)
% mirkeystrength   - Key strengths (probability of key candidates)
% mirkey           - Best keys and modes (in the 12 tone system)
% mirkeysom        - Visualizes key strengths with self-organizing map
% mirmode          - General estimation of mode (major/minor)
% mirtonalcentroid - Tonal centroid (using circles of fifths and thirds)
% mirhcdf          - Harmonic Change Detection Function
%
%RHYTHM
% mirtempo         - Tempo (in beats per minute)
% mirfluctuation   - Fluctuation strength (periodicities in each channel)
%
%TIMBRE
% mirmfcc          - Mel-frequency cepstrum coefficients
%                       (numerical description of the spectrum envelope)
% mirinharmonicity - Inharmonicity (partials non-multiple of fundamental)
% mirroughness     - Roughness (sensory dissonance)
% mirregularity    - Spectrum irregularity (amplitude variability of 
%                        successive peaks)
%
%ONSETS
% mironsets        - Note onset positions
% mirattacks       - Starting position of note attacks
% mirattacktime    - Duration of note attacks
% mirattackslope   - Average slope of note attacks
%
%ANALYSIS OF TEMPORAL FEATURES
% mirflux          - Flux, i.e., distance between successive frames
%
%ANALYSIS OF CURVES
% mirpeaks         - Peaks
% mirhisto         - Histogram
% mirentropy       - Entropy
% mirzerocross     - Sign-changes ratio
%
%ANALYSIS OF DISTRIBUTIONS (either spectra, envelopes, or histograms)
% mircentroid      - Centroid (center of gravity)
% mirspread        - Spread (non-concentration)
% mirskewness      - Skewness (lack of symmetry)
% mirkurtosis      - Kurtosis (peakiness)
% mirflatness      - Flatness
%
%SIMILARITY
% mirsimatrix      - Similarity matrix
% mirnovelty       - Novelty score
% mirdist          - Distance between audio files
% mirquery         - Query by example
%
%CLUSTERING
% mirclassify      - Classifies audio sequences
% mircluster       - Clusters segments in audio sequences
% mirhmm           - Sequence of states of a Hidden Markov Model
%
%MATLAB FUNCTIONS generalized to the MIRtoolbox data
% +                - Superposes audio files
% *                - Combines autocor, cepstrum curves
% corrcoef         - Computes correlation between curves
%
%DIVERSE
% mirchunklim      - Get or set the chunk size threshold
% mirwaitbar       - Toggles on/off the display of progress bars
% mirverbose       - Toggles on/off the display of ongoing operations


% mirpulseclarity  - Pulse clarity (strength of beats returned by mirtempo)
