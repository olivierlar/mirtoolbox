function t = mirtemporal(orig,varargin)
%   t = mirtemporal(x) creates a temporal object from signal x.
%   Optional arguments:
%       mirtemporal(...,'Center') centers the signal x.

if nargin > 0 && isa(orig,'mirtemporal')
    t.centered = orig.centered;
    t.nbits = orig.nbits;
else
    t.centered = 0;
    t.nbits = {};
end
t = class(t,'mirtemporal',mirdata(orig));
if nargin == 0 || not(isa(orig,'mirtemporal'))
    t = set(t,'Title','Temporal signal','Abs','time (s)','Ord','amplitude');
end
if nargin>1
    for i = 1:nargin-1
        if strcmp(varargin{i},'Center')
            d = get(t,'Data');
            for h = 1:length(d)
                for k = 1:length(d{h})
                    d{h}{k} = center(d{k});
                end
            end
            t = set(t,'Data',d);
        end
    end
    t = set(t,varargin{:});
end