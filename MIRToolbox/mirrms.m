function varargout = mirrms(x,varargin)
%   e = mirrms(x) calculates the root mean square energy.
%   Optional arguments:
%       mirrms(...,'Frame') computes the temporal evolution of the energy.

        notchunking.type = 'Boolean';
        notchunking.when = 'After';
        notchunking.default = 1;
    option.notchunking = notchunking;
    
        median.key = 'Median';
        median.when = 'Both';
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
if ~isamir(x,'mirscalar')
    if option.warning && ~isamir(x,'miraudio')
        warning(['Do you really intend to apply MIRRMS on a ',class(x),'?']);
    end
    d = get(x,'Data');
    v = mircompute(@algo,d,option.median);
    x = mirscalar(x,'Data',v,'Title','RMS energy');
end
if isstruct(postoption) && isfield(postoption,'notchunking') ...
        && postoption.notchunking
    if ~postoption.median
        x = after(x);
    end
elseif option.median
    mirerror('MIRRMS','''Median'' option should not be used in chunk decomposition mode. Results will not be reliable. Try mirchunklim(Inf)');
end


function e = algo(d,option_median)
nc = size(d,2);
nch = size(d,3);
e = zeros(1,nc,nch);
if option_median
    for i = 1:nch
        for j = 1:nc
            e(1,j,i) = sqrt(median(d(:,j,i).^2));
        end
    end
else
    for i = 1:nch
        for j = 1:nc
            e(1,j,i) = d(:,j,i)'*d(:,j,i);
        end
    end
end


function x = after(x)
v = mircompute(@afternorm,get(x,'Data'),get(x,'Length'));
x = set(x,'Data',v);

    
function d = afternorm(d,l)
d = sqrt(d/l);