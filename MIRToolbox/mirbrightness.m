function varargout = mirbrightness(x,varargin)
%   b = mirbrightness(s) calculates the spectral brightness, i.e. the amount
%       of spectral energy corresponding to frequencies higher than a given
%       cut-off threshold.
%   Optional arguments:
%   b = mirbrightness(s,'CutOff',f) specifies the frequency cut-off 
%       threshold in Hz.
%           Default value: f = 1500 Hz.
%
% Typical values for the frequency cut-off threshold:
%       3000 Hz in Juslin 2000, p. 1802.
%       1000 Hz and 500 Hz in Laukka, Juslin and Bresin 2005.
%
%   Juslin, P. N. (2000). Cue utilization in communication of emotion in 
% music performance: relating performance to perception. Journal of 
% Experimental Psychology: Human Perception and Performance, 26(6), 1797?813.
%   Laukka, P., Juslin, P. N., and Bresin, R. (2005). A dimensional approach 
% to vocal expression of emotion. Cognition and Emotion, 19, 633?653.


        cutoff.key = 'CutOff';
        cutoff.type = 'Integer';
        cutoff.default = 1500;
    option.cutoff = cutoff;
    
        minrms.key = 'MinRMS';
        minrms.when = 'After';
        minrms.type = 'Numerical';
        minrms.default = .005;
    option.minrms = minrms;
    
specif.option = option;
specif.defaultframelength = .05;
specif.defaultframehop = .5;


varargout = mirfunction(@mirbrightness,x,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if not(isamir(x,'mirspectrum'))
    x = mirspectrum(x);
end
type = 'mirscalar';


function b = main(s,option,postoption)
if iscell(s)
    s = s{1};
end
m = get(s,'Magnitude');
f = get(s,'Frequency');
w = warning('query','MATLAB:divideByZero');
warning('off','MATLAB:divideByZero');
v = mircompute(@algo,m,f,option.cutoff);
warning(w.state,'MATLAB:divideByZero');
b = mirscalar(s,'Data',v,'Title','Brightness');
if isstruct(postoption)
    b = after(s,b,postoption.minrms);
end


function v = algo(m,f,k)
if not(any(max(f)>k))
    warning('WARNING in MIRBRIGHTNESS: Frequency range of the spectrum too low for the estimation of brightness.');
end
sm = sum(m);
v = sum(m(f(:,1,1) > k,:,:)) ./ sm;


function b = after(s,b,minrms)
v = get(b,'Data');
if minrms
    r = mirrms(s,'Warning',0);
    v = mircompute(@trimrms,v,get(r,'Data'),minrms);
end
b = set(b,'Data',v);

    
function d = trimrms(d,r,minrms)
r = r/max(max(r));
pos = find(r<minrms);
for i = 1:length(pos)
    d(pos(i)) = NaN;
end
d = {d};