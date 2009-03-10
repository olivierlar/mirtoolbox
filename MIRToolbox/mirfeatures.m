function r = mirfeatures(x,varargin)
%   f = mirfeatures(x) computes a large set of features from one or several
%       audio files. x can be either the name of an audio file, or the
%       'Folder' keyword.
%   mirfeatures(...,'Stat') returns the statistics of the features instead
%       of the complete features themselves.
%   mirfeatures(...,'Segment') segments the audio sequence and analyzes
%       each segment separately.


%(not available yet)
%   mirfeatures(...,'Filterbank',nc) computes the analysis on each channel
%       of a filterbank decomposition.
%       Default value: nc = 5
%   mirfeatures(...,'Frame',...) 
%   mirfeatures(...,'Normal')
%   mirfeatures(...,'Sampling',s)



%RESAMPLING (22050)
%MIRAUDIO OPTIONS (Extract)

% CHECK BEHAVIOR WITH SHORT FILES

[stat,nchan,segm,feat] = scanargin(varargin);

if isa(x,'miraudio') || isa(x,'mirdesign')
    a = miraudio(x,'Normal'); % normalize with respect to RMS energy 
                              % in order to consider timbre independently of
                             % energy
else
    a = miraudio('Design','Normal');
    r = mirstruct;
end

if not(isempty(segm))
    a = mirsegment(a,segm);
end


% DYNAMICS
% --------

r.dynamics.rms = mirrms(a,'Frame');
% Perceived dynamics: spectral slope?

% RHYTHM
% ------

r.rhythm.tmp.fluctuation = mirfluctuation(a,'Summary');
r.rhythm.fluctuation.peak = mirpeaks(r.rhythm.tmp.fluctuation,'Total',1);%only one?
r.rhythm.fluctuation.centroid = mircentroid(r.rhythm.tmp.fluctuation);

r.rhythm.tmp.onsets = mironsets(a);

%r.rhythm.eventdensity = ...

r.rhythm.tempo = mirtempo(r.rhythm.tmp.onsets,'Frame');
%r.rhythm.pulseclarity = mirpulseclarity(r.tmp.onsets,'Frame');
    % Should use the second output of mirtempo.

r.rhythm.tmp.attacks = mironsets(r.rhythm.tmp.onsets,'Attacks');
r.rhythm.attack.time = mirattacktime(r.rhythm.tmp.attacks);
r.rhythm.attack.slope = mirattackslope(r.rhythm.tmp.attacks);

% TIMBRE
% ------

f = mirframe(a,.05,.5);
if 1 %max(strcmpi('centroid',feat)) || max(strcmpi('mfcc',feat))
    r.timbre.tmp.s = mirspectrum(f);
end
r.tmp.pitch = mirpitch(a,'Frame',.05,.5);

r.timbre.zerocross = mirzerocross(f);
if 1 %max(strcmpi('centroid',feat))
    r.timbre.centroid = mircentroid(r.timbre.tmp.s);
end
r.timbre.brightness = mirbrightness(r.timbre.tmp.s);
r.timbre.spread = mirspread(r.timbre.tmp.s);
r.timbre.skewness = mirskewness(r.timbre.tmp.s);
r.timbre.kurtosis = mirkurtosis(r.timbre.tmp.s);
r.timbre.rolloff95 = mirrolloff(r.timbre.tmp.s,95);
r.timbre.rolloff85 = mirrolloff(r.timbre.tmp.s,85);
r.timbre.spectentropy = mirentropy(r.timbre.tmp.s);
r.timbre.flatness = mirflatness(r.timbre.tmp.s);

r.timbre.roughness = mirroughness(r.timbre.tmp.s);
r.timbre.irregularity = mirregularity(r.timbre.tmp.s);
r.timbre.inharmonicity = mirinharmonicity(r.timbre.tmp.s,'f0',r.tmp.pitch);

if 1 %max(strcmpi('mfcc',feat))
    r.timbre.tmp.mfcc = mirmfcc(r.timbre.tmp.s);
    r.timbre.mfcc = mirmfcc(r.timbre.tmp.mfcc); %Direct assignment does not work yet...
    r.timbre.tmp.dmfcc = mirmfcc(r.timbre.tmp.mfcc,'Delta');
    r.timbre.dmfcc = mirmfcc(r.timbre.tmp.dmfcc); %Direct assignment does not work yet...
    r.timbre.ddmfcc = mirmfcc(r.timbre.tmp.dmfcc,'Delta');
end

r.timbre.lowenergy = mirlowenergy(f);
r.timbre.spectralflux = mirflux(f);

% PITCH
% -----

r.pitch.salient = mirpitch(r.tmp.pitch);%Direct assignment does not work yet...
r.tmp.chromagram = mirchromagram(a,'Frame','Wrap',0,'Pitch',0);
r.pitch.chromagram.peak=mirpeaks(r.tmp.chromagram,'Total',1);
r.pitch.chromagram.centroid=mircentroid(r.tmp.chromagram);

% TONALITY/HARMONY
% ----------------

r.tmp.keystrengths = mirkeystrength(r.tmp.chromagram);
[k1 ks]=mirkey(r.tmp.keystrengths,'Total',1);
r.tonal.keyclarity = ks;
r.tonal.mode = mirmode(r.tmp.keystrengths);
r.tonal.hcdf = mirhcdf(r.tmp.chromagram);

if stat
    r = mirstat(r);
    % SHOULD COMPUTE STAT OF CURVES FROM FRAMED_DECOMPOSED HIGH FEATURES
end
    
if not(isa(x,'miraudio')) && not(isa(x,'mirdesign'))
    r = mireval(r,x);
end


function [stat,nchan,segm,feat] = scanargin(v)
stat = 0;
nchan = 1;
segm = [];
feat = {};
i = 1;
while i <= length(v)
    arg = v{i};
    if ischar(arg) && strcmpi(arg,'Filterbank')
        i = i+1;
        if i <= length(v)
            nchan = v{i};
        else
            nchan = 10;
        end
    elseif ischar(arg) && strcmpi(arg,'Stat')
        i = i+1;
        if i <= length(v)
            stat = v{i};
        else
            stat = 1;
        end
    elseif ischar(arg) && strcmpi(arg,'Segment')
        i = i+1;
        if i <= length(v)
            segm = v{i};
        else
            segm = 1;
        end
    else
        feat{end+1} = arg;
    end    
    i = i+1;
end