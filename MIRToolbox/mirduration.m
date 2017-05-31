function varargout = mirduration(orig,varargin)
%   a = mirduration(x) estimates the duration of each note.
%   Optional arguments:
%   mirduration(...,'Contrast',c) specifies the 'Contrast' parameter
%       used in mironsets for event detection through peak picking.
%       Same default value as in mironsets.
%   mirduration(...,'Single') only selects one attack and decay phase in
%       the ignal (or in each segment).
    
        cthr.key = 'Contrast';
        cthr.type = 'Integer';
        cthr.default = .3;
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
    
specif.option = option;

varargout = mirfunction(@mirduration,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
o = mironsets(x,'Attack','Decay','Contrast',option.cthr,'Single',option.single,...
                 'Log',option.log,'MinLog',option.minlog,...
                'Presilence',option.presilence,'PostSilence',option.postsilence);
type = mirtype(x);


function du = main(o,option,postoption)
if iscell(o)
    o = o{1};
end
pa = get(o,'AttackPos');
pr = get(o,'DecayPos');
on = get(o,'OnsetPos');
off = get(o,'OffsetPos');
d = get(o,'Data');
t = get(o,'Pos');
[sl fp] = mircompute(@algo,pa,pr,d,t,on,off);
%fp = mircompute(@frampose,pru);
du = mirscalar(o,'Data',sl,'FramePos',fp,'Title','Duration','Unit','s.');
du = {du,o};


function fp = frampose(pr)
if isempty(pr)
    fp = [];
    return
end
pr = pr{1};
fp = pr;


function [du fp] = algo(pa,pr,d,t,on,off)
if isempty(pa)
    du = [];
    fp = [];
    return
end
pa = pa{1};
pr = pr{1};
on = on{1};
off = off{1};
du = zeros(1,length(pa));
fp = zeros(2,length(pa));
for i = 1:length(pa)
    [mv mp] = max(d(pa(i):pr(i)));
    mp = pa(i) + mp - 1;
    f1 = find(d(mp:-1:on(i)) < mv * .4,1);
    if isempty(f1)
        t1 = t(pa(i));
    else
        t1 = t(mp - f1);
    end
    f2 = find(d(mp:off(i)) < mv * .4,1);
    if isempty(f2)
        t2 = t(pr(i));
    else
        t2 = t(mp + f2);
    end
    du(i) = t2 - t1;
    fp(:,i) = [t1;t2];
end