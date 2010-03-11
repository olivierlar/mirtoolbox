function varargout = mirplay(a,varargin)
% mirplay method for mirmidi objects.

specif.option = struct;

specif.eachchunk = 'Normal';

varargout = mirfunction(@mirplay,a,varargin,nargout,specif,@init,@main);
if nargout == 0
    varargout = {};
end


function [x type] = init(x,option)
type = '';


function noargout = main(a,option,postoption)
if iscell(a)
    a = a{1};
end
d = get(a,'Data');
n = get(a,'Name');
for k = 1:length(d)
    display(['Playing analysis of file: ' n{k}])   
    playmidi(d{k});
end
noargout = {};