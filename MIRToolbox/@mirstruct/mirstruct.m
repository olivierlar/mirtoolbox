function d = mirstruct(varargin)

d.fields = {};
d.data = {};
d.tmp = struct;
d = class(d,'mirstruct');
d = set(d,varargin{:});