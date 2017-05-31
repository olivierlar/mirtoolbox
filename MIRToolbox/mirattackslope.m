function varargout = mirattackslope(orig,varargin)
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
    
        attack.key = {'Attack','Attacks'};
        attack.type = 'String';
        attack.choice = {'Derivate','Effort'};
        attack.default = 'Derivate';
    option.attack = attack;

        normal.key = 'Normal';
        normal.type = 'String';
        normal.choice = {0,1,'AcrossSegments'};
        normal.default = 'AcrossSegments';
    option.normal = normal;
    
        envmeth.type = 'String';
        envmeth.choice = {'Filter','Spectro'};
        envmeth.default = 'Spectro';
    option.envmeth = envmeth;

        ds.key = {'Down','PostDecim'};
        ds.type = 'Integer';
        if isamir(orig,'mirenvelope')
            ds.default = 1;
        else
            ds.default = NaN;
        end
    option.ds = ds;
    
specif.option = option;

varargout = mirfunction(@mirattackslope,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
if isnan(option.ds)
    if strcmpi(option.envmeth,'Spectro')
        option.ds = 0;
    else
        option.ds = 16;
    end
end
o = mironsets(x,option.envmeth,'Attack',option.attack,'Down',option.ds,...
                'Contrast',option.cthr,'Single',option.single,...
                'Log',option.log,'MinLog',option.minlog,...
                 option.envmeth,...
                'Presilence',option.presilence,'PostSilence',option.postsilence,...
                'Normal',option.normal);
type = mirtype(x);


function sl = main(o,option,postoption)
if iscell(o)
    o = o{1};
end
ap = get(o,'AttackPos');
op = get(o,'OnsetPos');
apu = get(o,'AttackPosUnit');
opu = get(o,'OnsetPosUnit');
sr = get(o,'Sampling');
d = get(o,'Data');
sl = mircompute(@algo,op,ap,opu,apu,d,option.meth,sr);
fp = mircompute(@frampose,opu,apu);
sl = mirscalar(o,'Data',sl,'FramePos',fp,'Title','Attack Slope');
sl = {sl,o};


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