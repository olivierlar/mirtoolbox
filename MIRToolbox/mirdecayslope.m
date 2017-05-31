function varargout = mirdecayslope(orig,varargin)
%   a = mirattackslope(x) estimates the average slope of each note decay.
%       Values are expressed in the same scale than the original signal,
%       but normalised by time in seconds.
%   Optional arguments:
%   a = mirdecayslope(x,m) specifies a method for slope computation.
%       Possible values:
%           m = 'Diff': ratio between the magnitude difference at the 
%               beginning and the ending of the attack period, and the
%               corresponding time difference.
%           m = 'Gauss': average of the slope, weighted by a gaussian
%               curve that emphasizes values at the middle of the attack
%               period. (similar to Peeters 2004).
%   mirdecayslope(...,'Contrast',c) specifies the 'Contrast' parameter
%       used in mironsets for event detection through peak picking.
%       Same default value as in mironsets.
%   mirdecayslope(...,'Single') only selects one attack phase in the
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

        log.key = {'LogOnset','LogCurve'};
        log.type = 'Boolean';
        log.default = 0;
    option.log = log;
    
        minlog.key = 'MinLog';
        minlog.type = 'Integer';
        minlog.default = 0;
    option.minlog = minlog;  
  
        envmeth.type = 'String';
        envmeth.choice = {'Filter','Spectro'};
        envmeth.default = 'Spectro';
    option.envmeth = envmeth;

        presilence.key = 'PreSilence';
        presilence.type = 'Boolean';
        presilence.default = 0;
    option.presilence = presilence;

        postsilence.key = 'PostSilence';
        postsilence.type = 'Boolean';
        postsilence.default = 0;
    option.postsilence = postsilence;

        normal.key = 'Normal';
        normal.type = 'String';
        normal.choice = {0,1,'AcrossSegments'};
        normal.default = 'AcrossSegments';
    option.normal = normal;

specif.option = option;

varargout = mirfunction(@mirdecayslope,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
o = mironsets(x,'Decay','Contrast',option.cthr,'Single',option.single,...
                 'Log',option.log,'MinLog',option.minlog,...
                 option.envmeth,...
                'Presilence',option.presilence,'PostSilence',option.postsilence,...
                'Normal',option.normal);
type = mirtype(x);


function rl = main(o,option,postoption)
if iscell(o)
    o = o{1};
end
rp = get(o,'DecayPos');
op = get(o,'OffsetPos');
rpu = get(o,'DecayPosUnit');
opu = get(o,'OffsetPosUnit');
sr = get(o,'Sampling');
d = get(o,'Data');
rl = mircompute(@algo,op,rp,opu,rpu,d,option.meth,sr);
fp = mircompute(@frampose,opu,rpu);
rl = mirscalar(o,'Data',rl,'FramePos',fp,'Title','Decay Slope');
rl = {rl,o};


function fp = frampose(op,ap)
if isempty(op)
    fp = [];
    return
end
op = sort(op{1});
ap = sort(ap{1});
fp = [op(:)';ap(:)'];


function sl = algo(po,pa,pou,pau,d,meth,sr)
if isempty(pa)
    sl = [];
    return
end
pa = sort(pa{1});
po = sort(po{1});
pau = sort(pau{1});
pou = sort(pou{1});
sl = zeros(1,length(pa));
for i = 1:length(pa)
    switch meth
        case 'Diff'
            sl(i) = (d(pa(i))-d(po(i)))/(pau(i)-pou(i));
        case 'Gauss'
            l = pa(i)-po(i);
            h = ceil(l/2);
            gauss = exp(-(1-h:l-h).^2/(l/4)^2);
            dat = diff(d(po(i):pa(i))).*gauss';
            sl(i) = mean(dat)*sr;
    end
end