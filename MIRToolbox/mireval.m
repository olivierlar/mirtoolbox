function v = mireval(d,file,single,export)
%   mireval(d,filename) applies the mirdesign object d to the audio file
%       named filename.
%   mireval(d,'Folder') applied the mirdesign object to each audio files in
%       the current directory.
%   mireval(d,'Folders') applied the mirdesign object recursively to the
%       subfolders.
%   Optional argument: mireval(...,'Single') only keeps the first
%       output when several output are returned for a given mirdesign
%       object.

% mireval performs the actual evaluation of the design flowchart.
%   If 'Folder' is used, the evaluation is carried out for each audio file
%       successively.
%   If d is a structure or a cell array, evaluate each component
%       separately.
%   The evaluation starts with a top-down traversal of the design flowchart
%       (evaleach).

if nargin<3
    single = [];
end
if nargin<4
    export = [];
end

% First, let's look at the content of the file(s): size, sampling rate,
% etc.
if not(ischar(file))
    error('ERROR IN EVAL: the second input should be a file name or ''Folder''')
end
w = [];    % Array containing the index positions of the starting and ending dates.
s = getsize(d);
if strcmpi(file,'Folder') || strcmpi(file,'Folders')
    [l w sr a] = evalfolder('',s,0,{},{},{},strcmpi(file,'Folders'));
    if l == 0
        disp('No sound file detected in this folder.')
    end
else
    l = 1;
    [d1,tp1,fp1,f1] = mirread([],file,0,0,0);
    if length(s)>1
        interval = s(1:2);
        if s(3)
            interval = round(interval*f1)+1;
        end
        if s(4) == 1
            interval = interval+round(d1/2);
        elseif s(4) == 2
            interval = interval+d1;
        end
        w = {min(max(interval,1),d1)};
    else
        w = {[1;d1]};
    end
    sr = {f1};
    a = {file};
end

if not(l)
    v = [];
    return
end

if isempty(export)
    y = cell(1,l);
end
order = 1:l;
if isa(d,'mirdesign') && isequal(get(d,'Method'),@mirplay)
    op = get(d,'Option');
    if isfield(op,'inc')
        if not(isnumeric(op.inc))
            op.inc = mirgetdata(op.inc);
        end
        [unused order] = sort(op.inc);
    elseif isfield(op,'dec')
        if not(isnumeric(op.inc))
            op.inc = mirgetdata(op.inc);
        end
        [unused order] = sort(op.dec,'descend');
    end
    if isfield(op,'every')
        order = order(1:op.every:end);
    end
    order = order(:)';
end

parallel = 0;
if mirparallel
    try
        matlabpool;
        parallel = 1;
        mirwaitbar(0)
        mirverbose(0)
    end
end

if parallel
    %   The evaluation is carried out for each audio file successively
    %       (or in parallel).
    parfor i = 1:l
        if l > 1
            fprintf('\n')
            display(['*** File # ',num2str(i),'/',num2str(l),': ',a{i}]);
        end
        tic
        yi = evalaudiofile(d,a{i},sr{i},w{i},{},0,i,single,'');
        toc
        y{i} = yi;
        if not(isempty(export))
            if strncmpi(export,'Separately',10)
                filename = a{i};
                filename(filename == '/') = '.';
                filename = [filename export(11:end)];
                if i == 1
                    mkdir('Backup');
                end
                mirexport(filename,yi);
            elseif i==1
                mirexport(export,yi);
            else
                mirexport(export,yi,'#add');
            end
        end
        clear yi
    end
else
    %   The evaluation is carried out for each audio file successively.
    isstat = isfield(d,'Stat');
    for i = 1:length(order)
        f = order(i);
        if l > 1
            fprintf('\n')
            display(['*** File # ',num2str(i),'/',num2str(l),': ',a{f}]);
        end
        tic
        yf = evalaudiofile(d,a{f},sr{f},w{f},{},0,f,single,'');
        toc
        y{f} = yf;
        if not(isempty(export))
            if strncmpi(export,'Separately',10)
                filename = a{f};
                filename(filename == '/') = '.';
                filename = ['Backup/' filename export(11:end)];
                if i == 1
                    mkdir('Backup');
                end
                mirexport(filename,yf);
            elseif i==1
                mirexport(export,yf);
            else
                mirexport(export,yf,'#add');
            end
        end
        clear yf
    end
end

%if isempty(export)
    v = combineaudiofile(a,isstat,y{:});
%else
%    v = [];
%end
    

function v = evalaudiofile(d,file,sampling,size,struc,istmp,index,single,name)
% Now let's perform the analysis (or analyses) on the different files.
%   If d is a structure or a cell array, evaluate each component
%       separately.
if isstruct(d)
    v = struct;
    if istmp
        struc.tmp = struct;
    end
    isstat = isfield(d,'Stat');
    if isstat
        d = rmfield(d,'Stat');
    end
    fields = fieldnames(d);
    for fi = 1:length(fields)
        fieldname = fields{fi};
        field = d.(fieldname);
        display(['*******',fieldname,'******']);
        if isstat
            if isa(field,'mirstruct')
                field = set(field,'Stat',1);
            elseif isa(field,'mirdesign')
                field = mirstat(field,'FileNames',0);
            else
                field.Stat = 1;
            end
        end
        res = evalaudiofile(field,file,sampling,size,struc,istmp,index,...
                                                        single,fieldname);
        if not(isempty(single)) && not(isequal(single,0)) && ...
                iscell(res) && isa(field,'mirdesign')
            res = res{1};
        end
        v.(fieldname) = res;
        if istmp
            struc.tmp.(fieldname) = res;
        end
        if fi == 1
            if isfield(res,'Class')
                v.Class = res.Class;
                v.(fieldname) = rmfield(res,'Class');
            end
        end
    end
    if isfield(v,'tmp')
        v = rmfield(v,'tmp');
    end
elseif iscell(d)
    l = length(d);
    v = cell(1,l);
    for j = 1:l
        v{j} = evalaudiofile(d{j},file,sampling,size,struc,istmp,index,single,[name,num2str(j)]);
    end
else
    d = set(d,'File',file,'Sampling',sampling,'Size',size,...
              'Eval',1,'Index',index,'Struct',struc);
    % For that particular file or this particular feature, let's begin the
    % actual evaluation process.
    v = evaleach(d,single,name);    
    % evaleach performs a top-down traversal of the design flowchart.
end


function c = combineaudiofile(filename,isstat,varargin) % Combine output from several audio files into one single
c = varargin{1};    % The (series of) input(s) related to the first audio file
if isempty(c)
    return
end
if isstruct(c)
    %fields = fieldnames(c);
    for j = 1:length(varargin)
        if j == 1
            fields = fieldnames(varargin{1});
        else
            fields = union(fields,fieldnames(varargin{j}));
        end
    end
    for i = 1:length(fields)
        field = fields{i};
        v = {};
        for j = 1:length(varargin)
            if isfield(varargin{j},field)
                v{j} = varargin{j}.(field);
            else
                v{j} = NaN;
            end
        end
        c.(field) = combineaudiofile('',isstat,v{:});
        if strcmp(field,'Class')
            c.Class = c.Class{1};
        end
    end
    if not(isempty(filename)) && isstat
        c.FileNames = filename;
    end
    return
end
if ischar(c)
    c = varargin;
    return
end
if (not(iscell(c)) && not(isa(c,'mirdata')))
    for j = 1:length(varargin)
        if j == 1
            lv = size(varargin{j},1);
        else
            lv = max(lv,size(varargin{j},1));
        end
    end
    c = NaN(lv,length(varargin));
    for i = 1:length(varargin)
        if not(isempty(varargin{i}))
            c(1:length(varargin{i}),i) = varargin{i};
        end
    end
    return
end
if (iscell(c) && not(isa(c{1},'mirdata')))
    for i = 1:length(c)
        for j = 1:nargin-2
            v{j} = varargin{j}{i};
        end
        c{i} = combineaudiofile(filename,isstat,v{:});
    end
    return
end
if not(iscell(c))
    c = {c};
end
nv = length(c); % The number of input variables for each different audio file
for j = 1:nv % Combine files for each different input variable
    v = varargin;
    for i = 1:length(varargin)
        if iscell(v{i})
            v{i} = v{i}{j};
        end
    end
    if not(isempty(v)) && not(isempty(v{1}))
        c{j} = combine(v{:});
    end
end


function s = getsize(d)
if isa(d,'mirstruct')
    d = get(d,'Data');
    if isempty(d)
        error('ERROR in MIREVAL: Your mirstruct object does not have any field (besides tmp).');
        s = 0;
    else
        s = getsize(d{1});
    end
elseif isstruct(d)
    fields = fieldnames(d);
    s = getsize(d.(fields{1}));
elseif iscell(d)
    s = getsize(d{1});
else
    s = get(d,'Size');  % Starting and ending dates in seconds.
end


function d2 = sortnames(d,d2,n)
if length(n) == 1
    d2(end+1) = d(1);
    return
end
first = zeros(1,length(n));
for i = 1:length(n)
    if isempty(n{i})
        first(i) = -Inf;
    else
        ni = n{i}{1};
        if ischar(ni)
            first(i) = ni-10058;
        else
            first(i) = ni;
        end
    end
end
[o i] = sort(first);
n = {n{i}};
d = d(i);
i = 0;
while i<length(n)
    i = i+1;
    if isempty(n{i})
        d2(end+1) = d(i);
    else
        dmp = (d(i));
        tmp = {n{i}(2:end)};
        while i+1<=length(n) && n{i+1}{1} == n{i}{1};
            i = i+1;
            dmp(end+1) = d(i);
            tmp{end+1} = n{i}(2:end);
        end
        d2 = sortnames(dmp,d2,tmp);
    end
end


function [l w sr a] = evalfolder(path,s,l,w,sr,a,folders)
if not(isempty(path))
    path = [path '/'];
end
dd = dir;
dn = {dd.name};
nn = cell(1,length(dn));  % Modified file names
for i = 1:length(dn)      % Each file name is considered
    j = 0;
    while j<length(dn{i})   % Each successive character is modified if necessary
        j = j+1;
        tmp = dn{i}(j) - '0';
        if tmp>=0 && tmp<=9
            while j+1<length(dn{i}) && dn{i}(j+1)>='0' && dn{i}(j+1)<='9'
                j = j+1;
                tmp = tmp*10 + (dn{i}(j)-'0');
            end
        else
            tmp = dn{i}(j);
        end
        nn{i}{end+1} = tmp;
    end
end
dd = sortnames(dd,[],nn);
for i = 1:length(dd);
    nf = dd(i).name;
    if folders && dd(i).isdir
        if not(strcmp(nf(1),'.'))
            cd(dd(i).name)
            [l w sr a] = evalfolder([path nf],s,l,w,sr,a,1);
            %l = l + l2;
            %w = [w w2];
            %sr = [sr sr2];
            %a = [a a2];
            cd ..
        end
    else
        [di,tpi,fpi,fi,bi,ni] = mirread([],nf,0,1,0);
        if not(isempty(ni))
            l = l+1;
            if not(isempty(s))
                interval = s(1:2);
                if s(3)
                    interval = round(interval*fi)+1;
                end
                if s(4) == 1
                    interval = interval+round(di/2);
                elseif s(4) == 2
                    interval = interval+di;
                end
                w{l} = min(max(interval,1),di);
            else
                w{l} = [1;di];
            end
            sr{l} = fi;
            a{l} = [path ni];
        end
    end
end