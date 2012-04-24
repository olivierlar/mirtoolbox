function varargout = mirattackleap(orig,varargin)
%   a = mirattackleap(x) estimates the leap of each note attack. 

        cthr.key = 'Contrast';
        cthr.type = 'Integer';
        cthr.default = NaN;
    option.cthr = cthr;

specif.option = option;

varargout = mirfunction(@mirattackleap,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
o = mironsets(x,'Attack','Contrast',option.cthr);
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
pa = sort(pa{1});
po = sort(po{1});
fp = [pa';po'];


function lp = algo(po,pa,d)
pa = sort(pa{1});
po = sort(po{1});
lp = zeros(1,length(pa));
for i = 1:length(pa)
    lp(i) = (d(po(i))-d(pa(i)));
end