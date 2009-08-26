function [y d2] = evaleach(d)
% Top-down traversal of the design flowchart, at the beginning of the
% evaluation phase.
% Called by mirfunction, mireval, mirframe and mirsegment.
% This is during that traversal that we check whether a chunk decomposition
% needs to be performed or not, and carry out that chunk decomposition.

CHUNKLIM = mirchunklim;
f = d.file;
fr = d.frame;
frnow = isfield(d.frame,'chunknow') && not(d.frame.chunknow);
sg = d.segment;
sr = d.sampling;
w = d.size;
lsz = w(2)-w(1)+1;    
len = lsz/sr;
if ischar(sg)
    error('ERROR in MIREVAL: mirsegment of design object accepts only array of numbers as second argument.');
end
if not(isempty(sg))
    over = find(sg > len,1);
    if not(isempty(over))
        sg = sg(1:over);
    end
end
a = d.argin;
ch = d.chunk;
specif = d.specif;
if isfield(specif,'combinechunk') && ...
        ischar(specif.combinechunk) && ...
        strcmpi(specif.combinechunk,'Average')
    specif.eachchunk = 'Normal';
end

if ischar(a)
    % The top-down traversal of the design flowchart now reaches the lowest
    % layer, i.e., audio file loading.
    % Now the actual evaluation will be carried out bottom-up.
    
    if isempty(ch)
        % No chunk decomposition
        y = miraudio(f,'Now',[w(:)' 0 0]);
    else
        % Chunk decomposition
        if isstruct(fr) && fr.eval
            % in a frame-decomposed context
            frs = fr.samples;
            y = miraudio(f,'Now',[frs(1,ch(1)),frs(2,ch(2)) 0]);
        else
            % in a non-frame-decomposed context            
            y = miraudio(f,'Now',[ch(1),ch(2) 0]);
        end
    end
    if not(isempty(d.postoption)) && d.postoption.mono
        y = miraudio(y,'Mono',1);
    end
    y = set(y,'AcrossChunks',get(d,'AcrossChunks'));
    if 0 %not(d.ascending)
        y = miraudio(y,'Reverse');
    end
    d2 = d;
    
elseif d.chunkdecomposed && isempty(d.tmpfile)
    % Already in a chunk decomposition process
    
    [y d2] = evalnow(d);  
    
elseif isempty(fr) || frnow || not(isempty(sg)) %% WHAT ABOUT CHANNELS?
    % No frame or segment decomposition in the design to evaluate
    % (Maybe implicit frame decomposition, though (frnow).)
    
    if not(isfield(specif,'eachchunk')) || d.nochunk
        chunks = [];
    elseif isempty(fr)
        if isempty(sg)
            meth = 'Chunk ';
            if lsz > CHUNKLIM
            % The required memory exceed the max memory threshold.
                nch = ceil(lsz/CHUNKLIM); 
            %%% TAKE INTO CONSIDERATION NUMBER OF CHANNELS; ETC... 
                chunks = max(1,lsz-CHUNKLIM*(nch:-1:1));
                chunks(2,:) = lsz-CHUNKLIM*(nch-1:-1:0)-1;
            else
                chunks = [];
            end
        else
            meth = 'Segment ';
            chunks = sg(1:end-1)*sr;
            chunks(2,:) = min( sg(2:end)*sr-1,lsz-1);
        end
    else
        meth = 'Chunk ';
        chunks = compute_frames(fr,sr,w,lsz,CHUNKLIM,d.overlap);
    end
    
    if not(isempty(chunks))
        % The chunk decomposition is performed.
        nch = size(chunks,2);
        d = callbeforechunk(d,d,w,lsz); % Some optional initialisation
        tmp = [];
        if mirwaitbar
            h = waitbar(0,['Computing ' func2str(d.method)]);
        else
            h = 0;
        end
        if not(isempty(d.tmpfile)) && d.tmpfile.fid == 0
            % When applicable, a new temporary file is created.
            d.tmpfile.fid = fopen('tmpfile.mirtoolbox','w');
        end
        tmpfile = [];
        if not(d.ascending)
            chunks = fliplr(chunks);
        end

        afterpostoption = d.postoption;
        method = d.method;
        if ischar(specif.eachchunk) && strcmpi(specif.eachchunk,'Normal')
            if not(isempty(d.postoption))
                pof = fieldnames(d.postoption);
                for o = 1:length(pof)
                    if isfield(specif.option.(pof{o}),'chunkcombine')
                        afterpostoption = rmfield(afterpostoption,pof{o});
                    else
                        d.postoption = rmfield(d.postoption,pof{o});
                    end
                end
            end
        else
            method = specif.eachchunk;
        end

        d2 = d;
        d2.method = method;
        for i = 1:size(chunks,2)
            disp([meth,num2str(i),'/',num2str(nch),'...'])
            d2 = set(d2,'Chunk',[chunks(1,i)+w(1) chunks(2,i)+w(1) (i == size(chunks,2))]);
            if not(ischar(specif.eachchunk) && ...
                   strcmpi(specif.eachchunk,'Normal'))
               if frnow
                   d2.postoption = 0;
               else
                    diffchunks = diff(chunks); % Usual chunk size
                    d2.postoption = max(diffchunks) -  diffchunks(i);
                        % Reduction of the current chunk size to be taken into
                        % consideration in mirspectrum, for instance, using
                        % zeropadding
               end
            end
            d2 = set(d2,'Tmp',tmp);
            d2.chunkdecomposed = 1;

            [ss d3] = evalnow(d2);
            % d2 is like d3 except that its argument is now evaluated.
            d3.postoption = d.postoption; % Pas joli joli
            d3.method = method;
            d2 = d3; % This new argument is transfered to d

            if isfield(specif,'combinechunk') && ...
                    ischar(specif.combinechunk) && ...
                    strcmpi(specif.combinechunk,'Average')
                % Measure total size for later averaging
                dss = get(ss,'Data');
                dss = mircompute(@multweight,dss,chunks(2,i)-chunks(1,i)+1);
                ss = set(ss,'Data',dss);
            end

            if isempty(ss)
                y = {};
            else
                if not(iscell(ss))
                    ss = {ss};
                end
                if isempty(sg)
                    tmp = get(ss{1},'Tmp');
                    if not(isempty(d2.tmpfile)) && d2.tmpfile.fid > 0
                        % If temporary file is used, chunk results are written
                        % in the file
                        if i < size(chunks,2)
                            ds = get(ss{1},'Data');
                            ps = get(ss{1},'Pos');
                            if 0 %(isfield(d.option,'zp') && d.option.zp == 2) ...
                                 %    || not(get(d,'Ascending'))
                                pos = chbeg*(size(ds{1}{1},3)+1);
                                if fseek(d2.tmpfile.fid,pos*8,'bof');
                                    nbl = fix(pos/2^16);
                                    rbl = rem(pos,2^16);
                                    for ibl = 1:nbl
                                        fwrite(d2.tmpfile.fid,zeros(2^16,1),'double');
                                    end
                                    if rbl
                                        fwrite(d2.tmpfile.fid,zeros(rbl,1),'double');
                                    end
                                end
                            end
                           % ftell(d2.tmpfile.fid)
                            count = fwrite(d2.tmpfile.fid,ds{1}{1},'double');
                            count = fwrite(d2.tmpfile.fid,ps{1}{1},'double');
                           % ftell(d2.tmpfile.fid)
                            clear ds ps
                        end
                        y = ss;
                    else
                        % Else, chunk results are directly combined in active
                        % memory
                        if i == 1
                            y = ss;
                        else
                            if isfield(specif,'combinechunk')
                                for z = 1:length(y)
                                    if ischar(specif.combinechunk)
                                        if strcmpi(specif.combinechunk,'Concat')
                                            y{z} = concatchunk(y{z},ss{z},d2.ascending);
                                        elseif strcmpi(specif.combinechunk,'Average')
                                            y{z} = sumchunk(y{z},ss{z});
                                        else
                                            error(['SYNTAX ERROR: ',...
                                                specif.combinechunk,...
                                        ' is not a known keyword for combinechunk.']);
                                        end
                                    else
                                        y{z} = specif.combinechunk(y{z},ss{z});
                                    end
                                end
                            else
                                y = {};
                            end
                        end
                    end
                else
                    if i == 1
                        y = ss;
                    else
                        for z = 1:length(y)
                            y{z} = combinesegment(y{z},ss{z});
                        end
                    end
                end
            end
            clear ss
            if h
                if not(d.ascending)
                    close(h)
                    h = waitbar((chunks(1,i)-chunks(1,end))/chunks(2,1),...
                        ['Computing ' func2str(d.method) ' (backward)']);
                else
                    waitbar((chunks(2,i)-chunks(1))/chunks(end),h)
                end
            end
        end
        if isfield(specif,'afterchunk')
            y{1} = specif.afterchunk(y{1},lsz,d.postoption);
        elseif isfield(specif,'combinechunk') && ...
                ischar(specif.combinechunk) && ...
                strcmpi(specif.combinechunk,'Average')
            y{1} = divideweightchunk(y{1},lsz);
        else
            if not(isempty(afterpostoption)) && isempty(d2.tmpfile)
                y{1} = d.method(y{1},[],afterpostoption);
            end
        end
        if not(isempty(d2.tmpfile))
            adr = ftell(d2.tmpfile.fid);
            fclose(d2.tmpfile.fid);
            ytmpfile.fid = fopen('tmpfile.mirtoolbox');
            fseek(ytmpfile.fid,adr,'bof');
            ytmpfile.data = y{1};
            ytmpfile.layer = 0;
            y{1} = set(y{1},'TmpFile',ytmpfile);
        end
        if h
            close(h)
        end
        drawnow

    else 
        % No chunk decomposition
        [y d2] = evalnow(d);
    end    
elseif d.nochunk
    [y d2] = evalnow(d);
else
    % No frame or segment decomposition in the design to be evaluated.
    [chunks fp] = compute_frames(fr,sr,w,lsz,CHUNKLIM,d.overlap);
    if size(chunks,2)>1
        % The chunk decomposition is performed.
        if mirwaitbar
            h = waitbar(0,['Computing ' func2str(d.method)]);
        else
            h = 0;
        end
        tmp = [];
        d = set(d,'Frames',fp);
        d2 = d;
        nch = size(chunks,2);
        for fri = 1:nch     % For each chunk...
            disp(['Chunk ',num2str(fri),'/',num2str(nch),'...'])
            d2 = set(d2,'Chunk',chunks(:,fri)');
            d2 = set(d2,'Tmp',tmp);
            %d2.postoption = [];
            [res d2] = evalnow(d2);
            if not(isempty(res))
                if iscell(res)
                    tmp = get(res{1},'Tmp');
                else
                    tmp = get(res,'Tmp');
                    res = {res};
                end
            end
            if fri == 1
                y = res;
            elseif isfield(specif,'combineframes')
                y = specif.combineframes(y{1},res{1});
            else
                y = combineframes(y,res);
            end
            if h
                waitbar(chend/nfr,h);
            end
        end
        if h
            close(h)
        end
    else
        % No chunk decomposition
        [y d2] = evalnow(d);
    end
end

 
if iscell(y)
    for i = 1:length(y)
        if not(isempty(y{i}))
            if iscell(y{i})
                for j = 1:length(y{i})
                    y{i}{j} = set(y{i}{j},'Tmp',[]);
                end
            else
                y{i} = set(y{i},'Tmp',[]);
            end
        end
    end
end


function [chunks fp] = compute_frames(fr,sr,w,lsz,CHUNKLIM,frov)
if strcmpi(fr.length.unit,'s')
    fl = ceil(fr.length.val*sr);
elseif strcmpi(fr.length.unit,'sp')
    fl = fr.length.val;
end
if strcmpi(fr.hop.unit,'/1')
    h = ceil(fr.hop.val*fl);
elseif strcmpi(fr.hop.unit,'%')
    h = ceil(fr.hop.val*fl*.01);
elseif strcmpi(fr.hop.unit,'s')
    h = ceil(fr.hop.val*sr);
elseif strcmpi(fr.hop.unit,'sp') || strcmpi(fr.hop.unit,'Hz')
    h = fr.hop.val;
end
if strcmpi(fr.hop.unit,'Hz')
    n = floor((lsz-fl)/sr*h)+1;   % Number of frames
else
    n = floor((lsz-fl)/h)+1;   % Number of frames
end
if n < 1
    %warning('WARNING IN MIRFRAME: Frame length longer than total sequence size. No frame decomposition.');
    fp = w;
else
    if strcmpi(fr.hop.unit,'Hz')
        st = floor(((1:n)-1)/h*sr)+w(1);
    else
        st = floor(((1:n)-1)*h)+w(1);
    end
    fp = [st; st+fl-1];
    fp(:,fp(2,:)>w(2)) = [];
end
fpsz = (fp(2,1)-fp(1,1)) * n;      % Total number of samples
if fpsz > CHUNKLIM
    % The required memory exceed the max memory threshold.
    nfr = size(fp,2);                     % Total number of frames
    frch = max(ceil(CHUNKLIM/fp(2,1)),2); % Number of frames per chunk
    frch = max(frch,frov*2);
    nch = ceil((nfr-frch)/(frch-frov))+1; % Number of chunks
    chbeg = (frch-frov)*(0:nch-1)+1;    % First frame in the chunk
    chend = (frch-frov)*(0:nch-1)+frch; % Last frame in the chunk
    chend = min(chend,nfr);
    if frov > 1
        chbeg = chend-frch+1;
    end
    chunks = [fp(1,chbeg) ; fp(2,chend)];
else
    chunks = [];
end


function old = combineframes(old,new)
if not(iscell(old))
    old = {old};
end
if not(iscell(new))
    new = {new};
end
for var = 1:length(new)
    ov = old{var};
    nv = new{var};
    if isa(ov,'mirscalar')
        ov = combinedata(ov,nv,'Data');
        ov = combinedata(ov,nv,'Mode');
    else
        if isa(ov,'mirtemporal')
            [ov omatch nmatch] = combinedata(ov,nv,'Time',[],[],@modiftime);
        else
            [ov omatch nmatch] = combinedata(ov,nv,'Pos',[],[]);
        end
        ov = combinedata(ov,nv,'Data',omatch,nmatch);
    end
    ov = combinedata(ov,nv,'FramePos');
    ov = combinedata(ov,nv,'PeakPos');
    ov = combinedata(ov,nv,'PeakVal');
    ov = combinedata(ov,nv,'PeakMode');
    old{var} = ov;
end


function [ov omatch nmatch] = combinedata(ov,nv,key,omatch,nmatch,modifdata)
odata = get(ov,key);
if isempty(odata) || isempty(odata{1})
    return
end
odata = odata{1};
if iscell(odata)
    if ischar(odata{1})
        return
    else
        odata = odata{1};
    end
end
ndata = get(nv,key);
ndata = ndata{1};
if iscell(ndata)
    ndata = ndata{1};
end
if nargin>3 
    if isempty(omatch)
        ol = size(odata,1);
        nl = size(ndata,1);
        unmatch = ol-nl;
        if unmatch>0
            [unused idx] = min(odata(1:1+unmatch,1,1)-ndata(1));
            omatch = idx:idx+nl-1;
            nmatch = 1:nl;
        elseif unmatch<0
            [unused idx] = min(ndata(1:1-unmatch,1,1)-odata(1));
            nmatch = idx:idx+ol-1;
            omatch = 1:ol;
        else
            nmatch = 1:nl;
            omatch = 1:ol;
        end       
    end
    odata(omatch,end+1:end+size(ndata,2),:,:) = ndata(nmatch,:,:,:); %4.D for keysom
else
    odata(:,end+1:end+size(ndata,2),:,:) = ndata;
end
ov = set(ov,key,{{odata}});  %{odata} for warped chromagram for instance....


function d = modiftime(d,p)
d = d + p;


function [y d] = evalnow(d)
% Go one step further in the top-down evaluation initialisation
argin = d.argin;
if not(iscell(argin))
    argin = {argin};
end
for i = 1:length(argin)
    a = argin{i};
    if not(d.ascending)
        a.ascending = 0;
    end
    if isa(a,'mirdata')
        % Input already computed
        tmpfile = get(a,'TmpFile');
        if not(isempty(tmpfile)) && tmpfile.fid > 0
            % The input can be read from the temporary file
            ch = get(d,'Chunk');
            a = tmpfile.data;
            a = set(a,'Tmp',get(d,'Tmp'),'TmpFile',tmpfile);
            channels = get(a,'Channels');
            channels = length(channels{1});
            if not(channels)
                channels = 1;
            end
            size = (ch(2)-ch(1)+1);
            current = ftell(tmpfile.fid);
            fseek(tmpfile.fid,current-size*(channels+1)*8,'bof');
            %ftell(tmpfile.fid)
            [data count] = fread(tmpfile.fid,[size,channels],'double');
            %count
            data = reshape(data,[size,1,channels]);
            [pos count] = fread(tmpfile.fid,size,'double');
            %count
           % ftell(tmpfile.fid)
            fseek(tmpfile.fid,current-size*(channels+1)*8,'bof');
            a = set(a,'Data',{{data}},'Pos',{{pos}});
            if ch(3)
                fclose(tmpfile.fid);
                delete('tmpfile.mirtoolbox');
            end
            argin{i} = a;
        end
    elseif isa(a,'mirdesign')
        if isempty(a.stored)
            % The design parameters are transfered to the previous component
            % in the design process
            a.size = d.size;
            a.chunk = d.chunk;
            a.file = d.file;
            a.eval = 1;
            a.tmp = d.tmp;
            a.sampling = d.sampling;
            if isstruct(d.frame) && isfield(d.frame,'samples') ...
                                 && not(isempty(d.frame.samples))
                a.chunkdecomposed = 1;
            else
                a.chunkdecomposed = d.chunkdecomposed;
            end
            if not(isempty(d.frame)) && ...
               not(strcmp(func2str(d.method),'mirframe'))
                a.frame = d.frame;
            end
            a.ready = 1;
            a.acrosschunks = d.acrosschunks;
            a.index = d.index;
            argin{i} = a;
        else
            % Variable already calculated
            tmp = get(d,'Struct');
            for j = 1:length(a.stored) % (if modified, modify also mirframe)
                stored = a.stored{j};
                if iscell(stored)
                    if length(stored)>1
                        tmp = tmp{stored{1},stored{2}};
                    else
                        tmp = tmp{stored{1}};
                    end
                else
                    tmp = getfield(tmp,stored);
                end
            end
            if iscell(tmp)
                tmp = tmp{1};
            end
            argin{i} = tmp;
        end
    end
end
if not(iscell(d.argin))
    argin = argin{1};
end
d.option.struct = get(d,'Struct');
if iscell(d.postoption)
    [y argin] = d.method(argin,d.option,d.postoption{:});
else
    [y argin] = d.method(argin,d.option,d.postoption);
end 
d = set(d,'Argin',argin);


function d0 = callbeforechunk(d0,d,w,lsz)
% If necessary, the chunk decomposition is performed a first time for
% initialisation purposes.
% Currently used only for miraudio(...,'Normal')
if not(ischar(d))
    specif = d.specif;
    CHUNKLIM = mirchunklim;
    nch = ceil(lsz/CHUNKLIM); 
    if isfield(specif,'beforechunk') ...
            && ((isfield(d.option,specif.beforechunk{2}) ...
                    && d.option.(specif.beforechunk{2})) ...
             || (isfield(d.postoption,specif.beforechunk{2}) ...
                    && d.postoption.(specif.beforechunk{2})) )
        if mirwaitbar
            h = waitbar(0,['Preparing ' func2str(d.method)]);
        else
            h = 0;
        end
        for i = 1:nch
            disp(['Chunk ',num2str(i),'/',num2str(nch),'...'])
            chbeg = CHUNKLIM*(i-1);
            chend = CHUNKLIM*i-1;
            d2 = set(d,'Size',d0.size,'File',d0.file,...
                       'Chunk',[chbeg+w(1) min(chend,lsz-1)+w(1)]);
            d2.method = specif.beforechunk{1};
            d2.postoption = {chend-lsz};
            d2.chunkdecomposed = 1;
            [tmp d] = evalnow(d2);
            d0 = set(d0,'AcrossChunks',tmp);
            if h
                waitbar(chend/lsz,h)
            end
        end
        if h
            close(h);
        end
        drawnow
    else
        d0 = callbeforechunk(d0,d.argin,w,lsz);
    end
end


function y = concatchunk(old,new,ascending)
do = get(old,'Data');
to = get(old,'Pos');
dn = get(new,'Data');
tn = get(new,'Pos');
if ascending
    y = set(old,'Data',{{[do{1}{1};dn{1}{1}]}},...
                'Pos',{{[to{1}{1};tn{1}{1}]}});
else
    y = set(old,'Data',{{[dn{1}{1};do{1}{1}]}},...
                'Pos',{{[tn{1}{1};to{1}{1}]}});
end


function y = combinesegment(old,new)

do = get(old,'Data');
to = get(old,'Pos');
fpo = get(old,'FramePos');
ppo = get(old,'PeakPos');
pppo = get(old,'PeakPrecisePos');
pvo = get(old,'PeakVal');
ppvo = get(old,'PeakPreciseVal');
pmo = get(old,'PeakMode');
apo = get(old,'AttackPos');
rpo = get(old,'ReleasePos');
tpo = get(old,'TrackPos');
tvo = get(old,'TrackVal');

dn = get(new,'Data');
tn = get(new,'Pos');
fpn = get(new,'FramePos');
ppn = get(new,'PeakPos');
pppn = get(new,'PeakPrecisePos');
pvn = get(new,'PeakVal');
ppvn = get(new,'PeakPreciseVal');
pmn = get(new,'PeakMode');
apn = get(new,'AttackPos');
rpn = get(new,'ReleasePos');
tpn = get(new,'TrackPos');
tvn = get(new,'TrackVal');

y = set(old,'Data',{{do{1}{:},dn{1}{:}}},...
            'FramePos',{{fpo{1}{:},fpn{1}{:}}}); 
        
if not(isempty(to))
    y = set(y,'Pos',{{to{1}{:},tn{1}{:}}}); 
end

if not(isempty(ppo))
    y = set(y,'PeakPos',{{ppo{1}{:},ppn{1}{:}}},...
                'PeakVal',{{pvo{1}{:},pvn{1}{:}}},...
                'PeakMode',{{pmo{1}{:},pmn{1}{:}}});
end

if not(isempty(pppn))
    y = set(y,'PeakPrecisePos',{[pppo{1},pppn{1}{1}]},...
                'PeakPreciseVal',{[ppvo{1},ppvn{1}{1}]});
end

if not(isempty(apn))
    y = set(y,'AttackPos',{[apo{1},apn{1}{1}]});
end

if not(isempty(rpn))
    y = set(y,'ReleasePos',{[rpo{1},rpn{1}{1}]});
end

if not(isempty(tpn))
    y = set(y,'TrackPos',{[tpo{1},tpn{1}{1}]});
end

if not(isempty(tvn))
    y = set(y,'TrackVal',{[tvo{1},tvn{1}{1}]});
end
 

function y = sumchunk(old,new,order)
do = mirgetdata(old);
dn = mirgetdata(new);
y = set(old,'ChunkData',do+dn);
        

function y = divideweightchunk(orig,length)
d = get(orig,'Data');
v = mircompute(@divideweight,d,length);
y = set(orig,'Data',v);

function e = multweight(d,length)
e = d*length;

function e = divideweight(d,length)
e = d/length;