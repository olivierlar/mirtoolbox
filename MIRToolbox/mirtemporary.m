function v = mirtemporary(s)
% mirtemporary(0) toggles off the progressive exportation of stats when
%   using the 'Folder' keyword
% mirtemporary(1) toggles back on the progressive exportation of stats.

persistent mir_temporary

if nargin
    mir_temporary = s;
else
    if isempty(mir_temporary)
        mir_temporary = 0;
    end
end

v = mir_temporary;