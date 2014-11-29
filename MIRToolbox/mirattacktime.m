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
    
        log.key = 'LogOnset';
        log.type = 'Boolean';
        log.default = 0;
    option.log = log;
    
        minlog.key = 'MinLog';
        minlog.type = 'Integer';
        minlog.default = 0;
    option.minlog = minlog;    

specif.option = option;

varargout = mirfunction(@mirattacktime,orig,varargin,nargout,specif,@init,@main);


function [o type] = init(x,option)
o = mironsets(x,'Attack','Contrast',option.cthr,'Single',option.single,...
                'Log',option.log,'MinLog',option.minlog,...
                'Filter','Normal','AcrossSegments');
type = mirtype(x);


function at = main(o,option,postoption)
if iscell(o)
    o = o{1};
end
po = get(o,'PeakPosUnit');
pa = get(o,'AttackPosUnit');
at = mircompute(@algo,po,pa,option.scale);
fp = mircompute(@frampose,pa,po);
at = mirscalar(o,'Data',at,'FramePos',fp,'Title','Attack Time');
at = {at,o};


function fp = frampose(pa,po)
if isempty(pa)
    fp = [];
    return
end
pa = sort(pa{1});
po = sort(po{1});
fp = [pa(:)';po(:)'];


function at = algo(po,pa,sc)
if isempty(pa)
    at = [];
    return
end
po = sort(po{1});
pa = sort(pa{1});
at = po-pa;
at = at';
if strcmpi(sc,'Log')
    at = log10(at);
end
at = at';