function varargout = mirdecreaseslope(varargin)
%   Same as mirreleaseslope

if nargout == 1
    varargout = {mirreleaseslope(varargin{:})};
else
    [rs,on] = mirreleaseslope(varargin{:});
    varargout = {rs,on};
end