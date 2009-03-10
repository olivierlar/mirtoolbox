function pb = mirwaitbar(s)

persistent mir_wait_bar

if nargin
    mir_wait_bar = s;
else
    if isempty(mir_wait_bar)
        mir_wait_bar = 1;
    end
end

pb = mir_wait_bar;