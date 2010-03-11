function p = mirparallel(s)
% mirparallel(1) toggles on parallel processing: when ?Folder? or ?Folders?
%   is used, several audio files can be analysed in parallel using several 
%   parallel Matlab sessions running on the different processors and/or 
%   processor cores of your computer.
%   (Requires MathWorks? Parallel Computing Toolbox.)
% mirparallel(0) toggles back off parallel processing.

persistent mir_parallel

if nargin
    mir_parallel = s;
else
    if isempty(mir_parallel)
        mir_parallel = 0;
    end
end

p = mir_parallel;