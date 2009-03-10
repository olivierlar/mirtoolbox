function s = mirscalar(orig,varargin)
%   s = mirscalar(x,n) creates a scalar object

if nargin == 0
    orig = [];
end
if iscell(orig)
    orig = orig{1};
end
if isa(orig,'mirscalar')
    s.mode = orig.mode;
    s.legend = orig.legend;
    s.parameter = orig.parameter;
else
    s.mode = [];
    s.legend = '';
    s.parameter = struct;
end
s = class(s,'mirscalar',mirdata(orig));
s = purgedata(s);
s = set(s,'Pos',{},'Abs','Temporal position of frames','Ord','Value',varargin{:});