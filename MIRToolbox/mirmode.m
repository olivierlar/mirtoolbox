function varargout = mirmode(x,varargin)
%   m = mirmode(a) estimates the mode. A value of 0 indicates a complete
%       incertainty, whereas a positive value indicates a dominance of
%       major mode and a negative value indicates a dominance of minor mode.
%   Optional arguments:
%   mirmode(a,s) specifies a strategy. 
%       Possible values for s: 'Sum', 'Best'(default)

        stra.type = 'String';
        stra.default = 'Best';
        stra.choice = {'Best','Sum','Major','SumBest'};
    option.stra = stra;
    
specif.option = option;
specif.defaultframelength = 1;
specif.defaultframehop = .5;

varargout = mirfunction(@mirmode,x,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if not(isamir(x,'mirkeystrength'))
    x = mirkeystrength(x);
end
type = 'mirscalar';


function o = main(s,option,postoption)
if iscell(s)
    s = s{1};
end
m = get(s,'Data');
v = mircompute(str2func(['algo' lower(option.stra)]),m);
b = mirscalar(s,'Data',v,'Title','Mode');
o = {b,s};


function v = algosum(m)
v = sum(abs(m(:,:,:,1) - m(:,:,:,2)));


function v = algobest(m)
v = max(m(:,:,:,1)) - max(m(:,:,:,2));


function v = algosumbest(m)
m = max(.5,m)-.5;
v = sum(m(:,:,:,1)) - sum(m(:,:,:,2));


function v = algomajor(m)
m = max(.5,m)-.5;
v = sum(m(:,:,1));