function varargout = mirduration(orig,varargin)
%   a = mirduration(x) estimates the duration of each note.
%   Optional arguments:
%   mirduration(...,'Contrast',c) specifies the 'Contrast' parameter
%       used in mironsets for event detection through peak picking.
%       Same default value as in mironsets.
%   mirduration(...,'Single') only selects one attack and release phase in
%       sthe ignal (or in each segment).
    
        cthr.key = 'Contrast';
        cthr.type = 'Integer';
        cthr.default = .3;
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

varargout = mirfunction(@mirduration,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
o = mironsets(x,'Attack','Release','Contrast',option.cthr,'Single',option.single,...
                 'Log',option.log,'MinLog',option.minlog,...
                 'Filter','Normal','AcrossSegments');
type = mirtype(x);


function du = main(o,option,postoption)
if iscell(o)
    o = o{1};
end
pa = get(o,'AttackPos');
pr = get(o,'ReleasePos');
d = get(o,'Data');
t = get(o,'Pos');
[sl fp] = mircompute(@algo,pa,pr,d,t);
%fp = mircompute(@frampose,pru);
du = mirscalar(o,'Data',sl,'FramePos',fp,'Title','Duration');
du = {du,o};


function fp = frampose(pr)
if isempty(pr)
    fp = [];
    return
end
pr = pr{1};
fp = pr;


function [du fp] = algo(pa,pr,d,t)
if isempty(pa)
    du = [];
    fp = [];
    return
end
pa = pa{1};
pr = pr{1};
du = zeros(1,length(pa));
fp = zeros(2,length(pa));
for i = 1:length(pa)
    [mv mp] = max(d(pa(i):pr(2,i)));
    mp = pa(i) + mp - 1;
    f1 = find(d(mp:-1:1) < mv * .4,1);
    if isempty(f1)
        t1 = t(pa(i));
    else
        t1 = t(mp - f1);
    end
    f2 = find(d(mp:pr(2,i)) < mv * .4,1);
    if isempty(f2)
        t2 = t(pr(2,i));
    else
        t2 = t(mp + f2);
    end
    du(i) = t2 - t1;
    fp(:,i) = [t1;t2];
end