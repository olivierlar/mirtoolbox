function s = mirstruct(varargin)

d = mirdesign('','',{},{},struct,{}); %,0);
s.fields = {};
s.data = {};
s.tmp = struct;
s.stat = 0;
s = class(s,'mirstruct',d);
s = set(s,varargin{:});