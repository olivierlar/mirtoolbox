function o = mirfunction(method,x,varg,nout,specif,init,main)

if isempty(x)
    o = {{},{},{}};
    return
end

if ischar(x) % Starting point of the design process
    design_init = 1;
    filename = x;
    orig = mirdesign(@miraudio,'Design',{varg},{},struct,'miraudio'); 
    % Implicitly, the audio file needs to be loaded first.
else
    design_init = 0;
    orig = x;
end

[orig during after] = miroptions(method,orig,specif,varg);

if isa(orig,'mirdesign')
    if not(get(orig,'Eval'))
        % Bottom-up construction of the general design
        
        if isstruct(during) && isfield(during,'frame') && ...
                isstruct(during.frame) && during.frame.auto
            % 'Frame' option: 
            % Automatic insertion of the mirframe step in the design
            orig = mirframe(orig,during.frame.length.val,...
                                 during.frame.length.unit,...
                                 during.frame.hop.val,...
                                 during.frame.hop.unit);   
        end
        
        % Automatic development of the implicit prerequisites,
        % with management of the data types throughout in the design process
        [orig type] = init(orig,during);
        
        if iscell(type)
            nout = length(type); %NUMBER OF OUTPUTS IS NOT MORE DEPENDENT ON THE CALL
        end
        
        o = mirdesign(method,orig,during,after,specif,type,nout);
        
        %if isstruct(during) && during.frame.auto
        %    % Now that the mirframe step has been integrated in the
        %    % flowchart, the 'Frame' option can be removed from the
        %    % flowchart's terminal node.
        %    opt = get(o,'Option');
        %    opt.frame = [];
        %    o = set(o,'Option',opt,'Frame',[]);
        %end
            
        if design_init && not(strcmpi(filename,'Design'))
            % When 'Design' keyword not used,
            % the function is immediately evaluated
            o = mireval(o,filename);
        else
            o = returndesign(o,nout);
        end
        if not(iscell(o))
            o = {o};
        end
        return
    else
        % Evaluation of the design.
        % First top-down initiation, then bottom-up process.
        
        if not(isempty(get(orig,'TmpFile'))) && get(orig,'ChunkDecomposed')
            %ch = get(orig,'Chunk');
            orig = evaleach(orig);
            if iscell(orig)
                orig = orig{1};
            end
            x = orig;
        else
            [orig x] = evaleach(orig);
        end
        
        %during.extract = [];
        if not(isequal(method,@nthoutput))
            if iscell(orig)
                orig = orig{1};
            end
            if isempty(get(orig,'Tmp'))
                orig = set(orig,'Tmp',get(x,'Tmp'));
            end
        end
    end
else
    design = 0;
    if iscell(orig)
        i = 0;
        while i<length(orig) && not(design)
            i = i+1;
            if isa(orig{i},'mirdesign')
                design = i;
            end
        end
    end
    if design
        % For function with multiple inputs
        if design == 1 && not(get(orig{1},'Eval'))
            % Progressive construction of the general design
            [orig type] = init(orig,during);
            o = mirdesign(method,orig,during,after,specif,type);
            o = set(o,'Size',get(orig{1},'Size'));
            o = returndesign(o,nout);
            return
        else
            % Progressive evaluation of the design
            for io = 1:length(orig)
                if isa(orig{io},'mirdesign')
                    o = evaleach(orig{io});
                    if iscell(o)
                        o = o{:};
                    end
                    orig{io} = o;
                end
            end
        end
    elseif not(isempty(init)) && not(isempty(during))
        if isstruct(during) && isfield(during,'frame') && ...
                isstruct(during.frame) && during.frame.auto
            orig = mirframe(orig,during.frame.length,during.frame.hop);        
        end
        orig = init(orig,during);
    end
end

% Actual computation
if not(iscell(orig) && not(ischar(orig{1}))) && ...
        not(isa(orig,'mirdesign') || isa(orig,'mirdata'))
    o = {orig};
    return
end
filenamearg = orig;
if iscell(filenamearg) && not(ischar(filenamearg{1}))
    filenamearg = filenamearg{1};
end
if iscell(filenamearg) && not(ischar(filenamearg{1}))
    filenamearg = filenamearg{1};
end
filename = get(filenamearg,'Name');
if not(isempty(during)) && mirverbose
    if length(filename) == 1
        disp(['Computing ',func2str(method),' related to ',filename{1},'...'])
    else
        disp(['Computing ',func2str(method),' for all audio files ...'])
    end
end
if not(iscell(orig) || isnumeric(x))
    orig = set(orig,'Index',get(x,'Index'));
end
o = main(orig,during,after);
if not(iscell(o) && length(o)>1) ...
        || (isa(x,'mirdesign') && get(x,'Eval')) ...
        || (iscell(x) && isa(x{1},'mirdesign') && get(x{1},'Eval'))
    o = {o x};
elseif (iscell(x) && isa(x{1},'mirdesign') && get(x{1},'Eval')) 
    o = {o x{1}};
end


function o = returndesign(i,nout)
o = cell(1,nout);
o{1} = i;
for k = 2:nout
    o{k} = nthoutput(i,k);
end