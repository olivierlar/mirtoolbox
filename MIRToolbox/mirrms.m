function varargout = mirrms(x,varargin)
%   e = mirrms(x) calculates the root mean square energy.
%   Optional arguments:
%       mirrms(...,'Frame') computes the temporal evolution of the energy.

        notchunking.type = 'Boolean';
        notchunking.when = 'After';
        notchunking.default = 1;
    option.notchunking = notchunking;
    
        median.key = 'Median';
        median.type = 'Boolean';
        median.default = 0;
    option.median = median;
        
        warning.key = 'Warning';
        warning.type = 'Boolean';
        warning.default = 1;
    option.warning = warning;

specif.option = option;

specif.defaultframelength = 0.05;
specif.defaultframehop = 0.5;

specif.combinechunk = 'Sum';

varargout = mirfunction(@mirrms,x,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
type = 'mirscalar';


function x = main(x,option,postoption)
if iscell(x)
    x = x{1};
end
if ~isempty(option) && option.median
    method = @median;
else
    method = @mean;
end
if ~isamir(x,'mirscalar')
    if option.warning && ~isamir(x,'miraudio')
        warning(['Do you really intend to apply MIRRMS on a ',class(x),'?']);
    end
    d = get(x,'Data');
    v = mircompute(@algo,d,method);
    x = mirscalar(x,'Data',v,'Title','RMS energy');
end
if isstruct(postoption) && isfield(postoption,'notchunking') ...
        && postoption.notchunking
    x = after(x);
end


function e = algo(d,method)
nc = size(d,2);
nch = size(d,3);
e = zeros(1,nc,nch);
for i = 1:nch
    for j = 1:nc
        dj = d(:,j,i).^2;
        e(1,j,i) = method(dj);
    end
end


function x = after(x)
v = mircompute(@afternorm,get(x,'Data'),get(x,'Length'));
x = set(x,'Data',v);

    
function d = afternorm(d,l)
d = sqrt(d);