function varargout = mirrolloff(x,varargin)
%   r = mirrolloff(s) calculates the spectral roll-off in Hz.
%   Optional arguments:
%   r = mirrolloff(s,'Threshold',p) specifies the energy threshold in 
%       percentage. (Default: .85)
%           p can be either a value between 0 and 1. But if p exceeds 1, it
%               is understood as a percentage, i.e. between 1 and 100.  
%           In other words, r is the frequency under which p percents
%               of the spectral energy is distributed.
%
% Typical values for the energy threshold:
%       85% in G. Tzanetakis, P. Cook. Musical genre classification of audio
%           signals. IEEE Tr. Speech and Audio Processing, 10(5),293-302, 2002.
%       95% in T. Pohle, E. Pampalk, G. Widmer. Evaluation of Frequently
%           Used Audio Features for Classification of Music Into Perceptual
%           Categories, ?

        p.key = 'Threshold';
        p.type = 'Integer';
        p.default = 85;
    option.p = p;

        minrms.key = 'MinRMS';
        minrms.when = 'After';
        minrms.type = 'Numerical';
        minrms.default = .005;
    option.minrms = minrms;

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
if option.p>1
    option.p = option.p/100;
end
v = mircompute(@algo,m,f,option.p);
r = mirscalar(s,'Data',v,'Title','Rolloff','Unit','Hz.');
if isstruct(postoption)
    r = after(s,r,postoption.minrms);
end


function v = algo(m,f,p)
cs = cumsum(m);          % accumulation of spectrum energy
thr = cs(end,:,:)*p;   % threshold corresponding to the rolloff point
v = zeros(1,size(cs,2),size(cs,3));
for l = 1:size(cs,3)
    for k = 1:size(cs,2)
        fthr = find(cs(:,k,l) >= thr(1,k,l),1); % find the location of the threshold
        if isempty(fthr)
            v(1,k,l) = NaN;
        else
            v(1,k,l) = f(fthr);
        end
    end
end


function r = after(s,r,minrms)
v = get(r,'Data');
if minrms
    rms = mirrms(s,'Warning',0);
    v = mircompute(@trimrms,v,get(rms,'Data'),minrms);
end
r = set(r,'Data',v);

    
function d = trimrms(d,r,minrms)
r = r/max(max(r));
pos = find(r<minrms);
for i = 1:length(pos)
    d(pos(i)) = NaN;
end
d = {d};