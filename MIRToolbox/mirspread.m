function varargout = mirspread(orig,varargin)
%   S = mirspread(x) calculates the spread of x, which can be either:
%       - a spectrum (spectral spread),
%       - an envelope (temporal spread), or
%       - any histogram.

        minrms.key = 'MinRMS';
        minrms.when = 'After';
        minrms.type = 'Numerical';
        minrms.default = .01;
    option.minrms = minrms;
    
specif.option = option;

varargout = mirfunction(@mirspread,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if not(isamir(x,'mirdata')) || isamir(x,'miraudio')
    x = mirspectrum(x);
end
type = 'mirscalar';


function S = main(x,option,postoption)
if iscell(x)
    x = x{1};
end
y = peaksegments(@spread,get(x,'Data'),...
                         get(x,'Pos'),...
                         get(mircentroid(x,'MaxEntropy',0),'Data'));
if isa(x,'mirspectrum')
    t = 'Spectral spread';
elseif isa(x,'mirenvelope')
    t = 'Temporal spread';
else
    t = ['Spread of ',get(x,'Title')];
end
S = mirscalar(x,'Data',y,'Title',t);
if isstruct(postoption)  && strcmpi(get(x,'Title'),'Spectrum') && ...
        isfield(postoption,'minrms') && postoption.minrms
    S = after(x,S,postoption.minrms);
end


function s = spread(d,p,c)
s = sqrt( sum((p-c).^2 .* (d/sum(d)) ) );


function S = after(x,S,minrms)
r = mirrms(x,'Warning',0);
v = mircompute(@trim,get(S,'Data'),get(r,'Data'),minrms);
S = set(S,'Data',v);

    
function d = trim(d,r,minrms)
r = r/max(max(r));
pos = find(r<minrms);
for i = 1:length(pos)
    d{pos(i)} = NaN;
end
d = {d};