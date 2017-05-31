function varargout = mirdecaytime(orig,varargin)
%   a = mirdecaytime(x) returns the duration (in s.) of each note decay. 
%   Optional arguments:
%   a = mirdecaytime(x,l) specifies whether to consider the duration in s.
%       (l='Lin') or the logarithm of that duration (l='Log') following the
%       approach proposed in Krimphoff et al. (1994).
%       Default value: l='Lin'.
%   mirdecaytime(...,'Single') only selects one attack phase in the signal
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

specif.option = option;

varargout = mirfunction(@mirdecaytime,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
o = mironsets(x,'Decay','Contrast',option.cthr,'Single',option.single,...
                'Log',option.log,'MinLog',option.minlog,...
                'Presilence',option.presilence,'PostSilence',option.postsilence,...
                'Normal','AcrossSegments');
type = mirtype(x);


function rt = main(o,option,postoption)
if iscell(o)
    o = o{1};
end
rp = get(o,'DecayPosUnit');
op = get(o,'OffsetPosUnit');
rt = mircompute(@algo,op,rp,option.scale);
fp = mircompute(@frampose,op,rp);
if strcmpi(option.scale,'Lin')
    unit = 's.';
else
    unit = '';
end
rt = mirscalar(o,'Data',rt,'FramePos',fp,'Title','Decay Time','Unit',unit);
rt = {rt,o};


function fp = frampose(op,rp)
if isempty(op)
    fp = [];
    return
end
op = sort(op{1});
rp = sort(rp{1});
fp = [op(:)';rp(:)'];


function rt = algo(op,rp,sc)
if isempty(rp)
    rt = [];
    return
end
op = sort(op{1});
rp = sort(rp{1});
rt = op-rp;
rt = rt';
if strcmpi(sc,'Log')
    rt = log10(rt);
end
rt = rt';