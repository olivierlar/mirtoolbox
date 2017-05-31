function varargout = mirattacktime(orig,varargin)
%   a = mirattacktime(x) returns the duration (in s.) of each note attack. 
%   Optional arguments:
%   a = mirattacktime(x,l) specifies whether to consider the duration in s.
%       (l='Lin') or the logarithm of that duration (l='Log') following the
%       approach proposed in Krimphoff et al. (1994).
%       Default value: l='Lin'.
%   mirattacktime(...,'Single') only selects one attack phase in the signal
%       (or in each segment).
%
% Krimphoff, J., McAdams, S. & Winsberg, S. (1994), Caractérisation du 
% timbre des sons complexes. II : Analyses acoustiques et quantification 
% psychophysique. Journal de Physique, 4(C5), 625-628.

        scale.type = 'String';
        scale.choice = {'Lin','Log'};
        scale.default = 'Lin';
    option.scale = scale;
    
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
    
        attack.key = {'Attack','Attacks'};
        attack.type = 'String';
        attack.choice = {'Derivate','Effort'};
        attack.default = 'Derivate';
    option.attack = attack;

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
    
        cutoff.key = 'CutOff';
        cutoff.type = 'Integer';
        cutoff.default = 37;
    option.cutoff = cutoff;
    
specif.option = option;

varargout = mirfunction(@mirattacktime,orig,varargin,nargout,specif,@init,@main);


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
                'Presilence',option.presilence,'PostSilence',option.postsilence,...
                'Normal','AcrossSegments','CutOff',option.cutoff);
type = mirtype(x);


function at = main(o,option,postoption)
if iscell(o)
    o = o{1};
end
ap = get(o,'AttackPosUnit');
op = get(o,'OnsetPosUnit');
at = mircompute(@algo,op,ap,option.scale);
fp = mircompute(@frampose,op,ap);
if strcmpi(option.scale,'Lin')
    unit = 's.';
else
    unit = '';
end
at = mirscalar(o,'Data',at,'FramePos',fp,'Title','Attack Time','Unit',unit);
at = {at,o};


function fp = frampose(op,ap)
if isempty(op)
    fp = [];
    return
end
op = sort(op{1});
ap = sort(ap{1});
fp = [op(:)';ap(:)'];


function at = algo(op,ap,sc)
if isempty(ap)
    at = [];
    return
end
op = sort(op{1});
ap = sort(ap{1});
at = ap-op;
at = at';
if strcmpi(sc,'Log')
    at = log10(at);
end
at = at';