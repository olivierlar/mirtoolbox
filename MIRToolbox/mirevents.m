function varargout = mirevents(varargin)
%   o = mirevents(x) shows a temporal curve where peaks relate to the
%       temporal position of events, and estimates those event time 
%       positions.
%   Optional arguments:
%       mirevents(...,f) selects the strategy for the computation of the
%           event detection function.
%           f = 'Envelope': Envelope of the audio signal. (Default choice).
%           With two methods for envelope extraction:
%               mirevents(...,'Spectro') (Default):
%                   mirevents(...,'SpectroFrame',fl,fh) species the frame
%                       length fl (in s.) and the hop factor fh (as a value
%                       between 0 and 1)
%                       Default values: fl = .1 s., fh = .1
%                    the frequency reassigment method can be specified:
%                    'Freq' (default), 'Mel', 'Bark' or 'Cents' (cf. mirspectrum).
%               mirevents(...,'Filter'):
%                   mirevents(...,'Filterbank',nc) specifies a preliminary
%                       filterbank decomposition into nc channels. If nc = 0,
%                       no decomposition is performed.
%                       Default value: 40.
%                   mirevents(...,'FilterbankType',ft) specifies the type of
%                       filterbank (see mirfilterbank).
%                       Default value: 'Gammatone';
%                   Options associated to the mirenvelope function can be
%                       passed here as well (see help mirenvelope):
%                      'FilterType','Tau','PreDecim'
%               mirevents(...,'Sum','no') does not sum back the channels at
%                   the end of the computation. The resulting detection curve
%                   remains therefore decomposed into several channels.
%               Options associated to the mirenvelope function can be
%                   passed here as well (see help mirenvelope):
%                   'HalfwaveCenter','Diff','HalfwaveDiff','Center',
%                   'Smooth', 'Sampling','Log','Power','Lambda',
%                  ,'PostDecim','UpSample'
%           f = 'SpectralFlux': Spectral flux of the audio signal.
%               Options associated to the mirflux function can be
%               passed here as well (see help mirflux):
%                   'Inc' (toggled on by default here),
%                   'Halfwave' (toggled on by default here),
%                   'Complex' (toggled off by default),
%                   'Median' (toggled on by default here)
%           f = 'Emerge': is an improved version of the 'SpectralFlux'
%               method that is able to detect more notes and in the same 
%               time ignore the spectral variation produced by vibrato.
%%%%
%   When the 'Emerge' method is used for academic research, please cite the 
%       following publication:
%   Lartillot, O., Cereghetti, D., Eliard, K., Trost, W. J., Rappaz, M.-A.,
%       Grandjean, D., "Estimating tempo and metrical features by tracking 
%       the whole metrical hierarchy", 3rd International Conference on 
%       Music & Emotion, Jyv?skyl?, 2013.
%%%%
%           f = 'Pitch ':computes a frame-decomposed autocorrelation function ,
%                of same default characteristics than those returned
%                by mirpitch, with however a range of frequencies set by 
%                the following options:
%                   'Min' (set by default to 30 Hz),
%                   'Max' (set by default to 1000 Hz),
%                and subsequently computes the novelty curve of the 
%                resulting similatrix matrix.
%               Option associated to the mirnovelty function can be
%               passed here as well (see help mirnovelty):
%                   'KernelSize' (set by default to 32 samples)
%       mirevents(...,'Detect',d) toggles on or off the event detection, 
%           which is based on the detection function.
%           (By default toggled on.)
%           Option associated to the mirpeaks function can be specified as
%               well:
%               'Contrast' with default value c = .01
%               'Threshold' with default value t = 0
%               'Single' detects only the highest peak.
%       mirevents(...,'Attack') (or 'Attacks') detects attack phases.
%       mirevents(...,'Decay') (or 'Decays') detects decay phases.
%       mirevents(...,'Frame',...) decomposes into frames, with default frame
%           length 3 seconds and hop factor .1
%   Preselected detection models:
%       mirevents(...,'Scheirer') corresponds to (Scheirer, 1998):
%           mirevents(...,'FilterBankType','Scheirer',...
%                         'FilterType','HalfHann','Sampling',200,...
%                         'HalfWaveDiff','Sum',0,'Detect',0)
%       mirevents(...,'Klapuri99') corresponds to most of (Klapuri, 1999).

varargout = {mironsets(varargin{:})};