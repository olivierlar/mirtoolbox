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

specif.option = option;

varargout = mirfunction(@mirattackleap,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
o = mironsets(x,'Attack','Contrast',option.cthr,'Single',option.single,...
                'Log',option.log,'MinLog',option.minlog,...
                'Filter','Normal','AcrossSegments');
type = mirtype(x);


function sl = main(o,option,postoption)
if iscell(o)
    o = o{1};
end
po = get(o,'PeakPos');
pa = get(o,'AttackPos');
pou = get(o,'PeakPosUnit');
pau = get(o,'AttackPosUnit');
d = get(o,'Data');
sl = mircompute(@algo,po,pa,d);
fp = mircompute(@frampose,pau,pou);
sl = mirscalar(o,'Data',sl,'FramePos',fp,'Title','Attack Leap');
sl = {sl,o};


function fp = frampose(pa,po)
if isempty(pa)
    fp = [];
    return
end
pa = sort(pa{1});
po = sort(po{1});
fp = [pa(:)';po(:)'];


function lp = algo(po,pa,d)
if isempty(pa)
    lp = [];
    return
end
pa = sort(pa{1});
po = sort(po{1});
lp = zeros(1,length(pa));
for i = 1:length(pa)
    lp(i) = (d(po(i))-d(pa(i)));
end