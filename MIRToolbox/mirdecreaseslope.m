function varargout = mirdecreaseslope(varargin)
%   Same as mirdecayslope

if nargout == 1
    varargout = {mirdecayslope(varargin{:})};
else
    [rs,on] = mirdecayslope(varargin{:});
    varargout = {rs,on};
end