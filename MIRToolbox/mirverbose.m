function v = mirverbose(s)

persistent mir_verbose

if nargin
    mir_verbose = s;
else
    if isempty(mir_verbose)
        mir_verbose = 1;
    end
end

v = mir_verbose;