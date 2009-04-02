function varargout = mirfluctuation(orig,varargin)
%   f = mirfluctuation(x) calculates the fluctuation strength, indicating the
%       rhythmic periodicities along the different channels.
%   Optional arguments:
%       mirfluctuation(...,'Summary') returns the summary of the fluctuation,
%           i.e., the summation along the critical bands.
%
% E. Pampalk, A. Rauber, D. Merkl, "Content-based Organization and 
% Visualization of Music Archives", 

        sum.key = 'Summary';
        sum.type = 'Boolean';
        sum.default = 0;
    option.sum = sum;

specif.option = option;
     
varargout = mirfunction(@mirfluctuation,orig,varargin,nargout,specif,@init,@main);


function [f type] = init(x,option)
if iscell(x)
    x = x{1};
end
if isamir(x,'miraudio') && not(isframed(x))
    x = mirframe(x,.023,.5);
end
m = mirspectrum(x,'Power','Terhardt','Bark','dB','Mask');
f = mirspectrum(m,'AlongBands','Max',10,'Window',0,...
                  'Resonance','Fluctuation','MinRes',.01);
if option.sum
    f = mirsummary(f);
end
type = 'mirspectrum';


function f = main(f,option,postoption)
f = set(f,'Title','Fluctuation');