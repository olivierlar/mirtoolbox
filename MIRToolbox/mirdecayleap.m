function varargout = mirdecayleap(orig,varargin)
%   a = mirdecayleap(x) estimates the leap of each note decay.
%       Values are expressed in the same scale than the original signal.
%   Optional arguments:
%   mirattackleap(...,'Contrast',c) specifies the 'Contrast' parameter
%       used in mironsets for event detection through peak picking.
%       Same default value as in mironsets.
%   mirattackleap(...,'Single') only selects one attack phase in the signal
%       (or in each segment).

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

varargout = mirfunction(@mirdecayleap,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
o = mironsets(x,'Decay','Contrast',option.cthr,'Single',option.single,...
                'Log',option.log,'MinLog',option.minlog,...
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
d = get(o,'Data');
rl = mircompute(@algo,op,rp,d);
fp = mircompute(@frampose,opu,rpu);
rl = mirscalar(o,'Data',rl,'FramePos',fp,'Title','Decay Leap');
rl = {rl,o};


function fp = frampose(po,pr)
if isempty(po) || isempty(pr)
    fp = [];
    return
end
po = sort(po{1});
pr = sort(pr{1});
fp = [po(:)';pr(:)'];


function lp = algo(po,pr,d)
if isempty(po) || isempty(pr)
    lp = [];
    return
end
po = sort(po{1});
pr = sort(pr{1});
lp = zeros(1,length(po));
for i = 1:length(po)
    lp(i) = (d(pr(i))-d(po(i)));
end