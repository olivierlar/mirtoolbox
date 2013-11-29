function varargout = mirdecreaseslope(orig,varargin)
%   a = mirattackslope(x) estimates the average slope of each note attack.
%       Values are expressed in the same scale than the original signal,
%       but normalised by time in seconds.
%   Optional arguments:
%   a = mirattackslope(x,m) specifies a method for slope computation.
%       Possible values:
%           m = 'Diff': ratio between the magnitude difference at the 
%               beginning and the ending of the attack period, and the
%               corresponding time difference.
%           m = 'Gauss': average of the slope, weighted by a gaussian
%               curve that emphasizes values at the middle of the attack
%               period. (similar to Peeters 2004).
%   mirattackslope(...,'Contrast',c) specifies the 'Contrast' parameter
%       used in mironsets for event detection through peak picking.
%       Same default value as in mironsets.
%   mirattackslope(...,'Single') only selects one attack phase in the
%       signal (or in each segment).
%
% Peeters. G. (2004). A large set of audio features for sound description
% (similarity and classification) in the CUIDADO project. version 1.0

        meth.type = 'String';
        meth.choice = {'Diff','Gauss'};
        meth.default = 'Diff';
    option.meth = meth;
    
        cthr.key = 'Contrast';
        cthr.type = 'Integer';
        cthr.default = NaN;
    option.cthr = cthr;
    
        single.key = 'Single';
        single.type = 'Boolean';
        single.default = 0;
    option.single = single;

        log.key = 'LogOnset';
        log.type = 'Boolean';
        log.default = 0;
    option.log = log;
    
        minlog.key = 'MinLog';
        minlog.type = 'Integer';
        minlog.default = 0;
    option.minlog = minlog;    

specif.option = option;

varargout = mirfunction(@mirdecreaseslope,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
o = mironsets(x,'Release','Contrast',option.cthr,'Single',option.single,...
                 'Log',option.log,'MinLog',option.minlog,...
                 'Filter','Normal','AcrossSegments');
type = mirtype(x);


function sl = main(o,option,postoption)
if iscell(o)
    o = o{1};
end
pr = get(o,'ReleasePos');
pru = get(o,'ReleasePosUnit');
sr = get(o,'Sampling');
d = get(o,'Data');
sl = mircompute(@algo,pr,pru,d,option.meth,sr);
fp = mircompute(@frampose,pru);
sl = mirscalar(o,'Data',sl,'FramePos',fp,'Title','Decrease Slope');
sl = {sl,o};


function fp = frampose(pr)
if isempty(pr)
    fp = [];
    return
end
pr = pr{1};
fp = pr;


function sl = algo(pr,pru,d,meth,sr)
if isempty(pr)
    sl = [];
    return
end
pr = pr{1};
pru = pru{1};
sl = zeros(1,length(pr));
for i = 1:size(pr,2)
    switch meth
        case 'Diff'
            sl(i) = (d(pr(1,i))-d(pr(2,i)))/(pru(2,i)-pru(1,i));
        case 'Gauss'
            l = pr(2,i)-pr(1,i);
            h = ceil(l/2);
            gauss = exp(-(1-h:l-h).^2/(l/4)^2);
            dat = -diff(d(pr(1,i):pr(2,i))).*gauss';
            sl(i) = mean(dat)*sr;
    end
end