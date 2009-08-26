function d = mirdesign(method,argin,option,postoption,specif,type,nout)

d.method = method;
d.argin = argin;
d.option = option;
d.postoption = postoption;
d.specif = specif;
d.type = type;
if ischar(argin)
    d.frame = {};
    d.segment = {};
    d.chunkdecomposed = 0;
    d.size = {};
    d.file = '';
    d.sampling = 0;
    d.nochunk = 0;
    if strcmp(func2str(method),'mirenvelope') && d.option.zp == 2
        %if ((isnan(d.option.zp) && strcmpi(option.filter,'IIR')) || ...
        %     (not(isnan(d.option.zp)) &&  d.option.zp)) 
        %    d.ascending = 1;
        %else
            d.ascending = 0;
        %end
    else
        d.ascending = 1;
    end
    d.overlap = 0;
else
    if iscell(argin)
        argin = argin{1};
    end
    if (strcmp(func2str(method),'mirspectrum') && d.option.alongbands) ...
        || (isfield(specif,'nochunk') && specif.nochunk)
        d.frame = [];
        if isfield(d.specif,'eachchunk')
            d.specif = rmfield(d.specif,'eachchunk');
            d.specif = rmfield(d.specif,'combinechunk');
        end
    else
        d.frame = argin.frame;
        if strcmp(func2str(method),'mirenvelope') ...
            && strcmpi(d.option.method,'Spectro') ...
            && not(isempty(d.frame))
                d.frame.chunknow = 0;
        end
    end
    d.segment = argin.segment;
    d.chunkdecomposed = argin.chunkdecomposed;
    d.size = argin.size;
    d.file = argin.file;
    d.sampling = argin.sampling;
    if (isfield(specif,'nochunk') && specif.nochunk) || not(isempty(argin.stored))
        d.nochunk = 1;
    else
        d.nochunk = argin.nochunk;
    end
    if strcmp(func2str(method),'mirenvelope')
        if d.option.zp == 2
            d.ascending = not(isempty(d.segment));
        else
            d.ascending = 1;
        end
    else
        d.ascending = argin.ascending;
    end
    d.overlap = argin.overlap;
end
d.chunk = [];
d.eval = 0;
d.tmp = [];
d.acrosschunks = []; % Data that can be accumulated among chunks during the beforechunk process.
d.ready = 0;
d.struct = [];
d.stored = [];
d.index = NaN;
if strcmp(func2str(method),'mirenvelope') && ...
            d.option.zp == 2 && isempty(d.segment)
    % Triggers the use of temporary file for the mirenvelope computation
    d.tmpfile.fid = 0;
else
    d.tmpfile = [];
end
if nargin < 7
    d.nout = 1;
else
    d.nout = nout;
end
d = class(d,'mirdesign');