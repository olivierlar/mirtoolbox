function varargout = mirattackleap(orig,varargin)
%   a = mirattackleap(x) estimates the leap of each note attack.
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

        log.key = 'LogOnset';
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
    
        attack.key = {'Attack','Attacks'};
        attack.type = 'String';
        attack.choice = {'Derivate','Effort'};
        attack.default = 'Derivate';
    option.attack = attack;

specif.option = option;

varargout = mirfunction(@mirattackleap,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
o = mironsets(x,'Attack',option.attack,'Contrast',option.cthr,'Single',option.single,...
                'Log',option.log,'MinLog',option.minlog,...
                'Presilence',option.presilence,'PostSilence',option.postsilence,...
                'Normal',option.normal);
type = mirtype(x);


function al = main(o,option,postoption)
if iscell(o)
    o = o{1};
end
ap = get(o,'AttackPos');
op = get(o,'OnsetPos');
apu = get(o,'AttackPosUnit');
opu = get(o,'OnsetPosUnit');
d = get(o,'Data');
al = mircompute(@algo,op,ap,d);
fp = mircompute(@frampose,opu,apu);
al = mirscalar(o,'Data',al,'FramePos',fp,'Title','Attack Leap');
al = {al,o};


function fp = frampose(po,pa)
if isempty(po)
    fp = [];
    return
end
po = sort(po{1});
pa = sort(pa{1});
fp = [po(:)';pa(:)'];


function lp = algo(po,pa,d)
if isempty(po)
    lp = [];
    return
end
po = sort(po{1});
pa = sort(pa{1});
lp = zeros(1,length(po));
for i = 1:length(po)
    lp(i) = (d(pa(i))-d(po(i)));
end