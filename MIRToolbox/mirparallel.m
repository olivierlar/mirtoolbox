function p = mirparallel(s)
% mirparallel(1) toggles on parallel processing.
%   When ?Folder? or ?Folders? is used, several audio files can be analysed
%   in parallel using several  parallel Matlab sessions running on the 
%   different processors and/or  processor cores of your computer.
%   (Requires MathWorks? Parallel Computing Toolbox.)
% mirparallel(0) toggles back off parallel processing.

persistent mir_parallel

if nargin
    if s
        try
            matlabpool size;
        catch
            mirerror('mirparallel','mirparallel does not work for Matlab 2013b and more recent. If you would like MIRtoolbox to be developed for parallel processing or for other purposes, and if you have any funding to suggest, please contact us.')
        end
    end
    mir_parallel = s;
else
    if isempty(mir_parallel)
        mir_parallel = 0;
    end
end

p = mir_parallel;