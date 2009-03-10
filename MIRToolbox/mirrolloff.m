function varargout = mirrolloff(x,varargin)
%   r = mirrolloff(s) calculates the spectral roll-off in Hz.
%   Optional arguments:
%   r = mirrolloff(s,'Threshold',p) specifies the energy threshold in 
%       percentage. (Default: .85)
%           In other words, r is the frequency under which 100*p percent
%               of the spectral energy is distributed.
%
% Typical values for the energy threshold:
%       85 in G. Tzanetakis, P. Cook. Musical genre classification of audio
%           signals. IEEE Tr. Speech and Audio Processing, 10(5),293-302, 2002.
%       95 in T. Pohle, E. Pampalk, G. Widmer. Evaluation of Frequently
%           Used Audio Features for Classification of Music Into Perceptual
%           Categories, ?

        p.key = 'Threshold';
        p.type = 'Integer';
        p.default = 85;
        p.position = 2;
    option.p = p;
    
specif.option = option;

varargout = mirfunction(@mirrolloff,x,varargin,nargout,specif,@init,@main);


function [s type] = init(x,option)
s = mirspectrum(x);
type = 'mirscalar';


function r = main(s,option,postoption)
if iscell(s)
    s = s{1};
end
m = get(s,'Magnitude');
f = get(s,'Frequency');
v = mircompute(@algo,m,f,option.p);
r = mirscalar(s,'Data',v,'Title','Rolloff','Unit','Hz.');


function v = algo(m,f,p)
cs = cumsum(m);          % accumulation of spectrum energy
thr = cs(end,:,:)*p/100;   % threshold corresponding to the rolloff point
v = zeros(1,size(cs,2),size(cs,3));
for l = 1:size(cs,3)
    for k = 1:size(cs,2)
        fthr = find(cs(:,k,l) >= thr(1,k,l)); % find the location of the threshold
        if isempty(fthr)
            v(1,k,l) = NaN;
        else
            v(1,k,l) = f(fthr(1));
        end
    end
end