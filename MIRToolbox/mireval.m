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

if not(ischar(file) || (iscell(file) && ischar(file{1})))
    error('ERROR IN MIREVAL: the second input should be a file name or ''Folder''')
end

if nargin<3
    single = [];
end
if nargin<4
    export = [];
end

% First, let's look at the content of the file(s): size, sampling rate,
% etc.
w = [];    % Array containing the index positions of the starting and ending dates.
s = getsize(d);
ch = 1;
if ~iscell(file) && (strcmpi(file,'Folder') || strcmpi(file,'Folders'))
    [l w sr lg a] = evalfolder('',s,0,[],[],[],{},strcmpi(file,'Folders'));
    if l == 0
        disp('No sound file detected in this folder.')
    end
elseif iscell(file) || (length(file)>3 && strcmpi(file(end-3:end),'.txt'))
    if iscell(file)
        a = file;
    else
        a = importdata(file);
    end
    l = length(a);
    w = zeros(2,l);
    sr = zeros(1,l);
    lg = zeros(1,l);
    for i = 1:l
        [di,tpi,fpi,fi,lg(i)] = mirread([],a{i},0,0,0);
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
            w(:,i) = min(max(interval,1),di);
        else
            w(:,i) = [1;di];
        end
        if getsampling(d)
            sr(i) = getsampling(d);
        else
            sr(i) = fi;
        end
    end
else
    l = 1;
    [d1,tp1,fp1,f1,lg,b,n,ch] = mirread([],file,0,0,0);
    if length(s)>1
        interval = s(1:2)';
        if s(3)
            interval = round(interval*f1)+1;
        end
        if s(4) == 1
            interval = interval+round(d1/2);
        elseif s(4) == 2
            interval = interval+d1;
        end
        if d1 < interval(2)
            warning('WARNING IN MIRAUDIO: The temporal region to be extracted exceeds the temporal extent of the whole audio file.'); 
        end
        w = min(max(interval,1),d1);
    else
        w = [1;d1];
    end
    if isa(d,'mirdesign') && getsampling(d)
        sr = getsampling(d);
    else
        sr = f1;
    end
    a = {file};
end

if not(l)
    v = [];
    return
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
        if mirparallel == 1 || mirparallel == Inf
            matlabpool;
        else
            matlabpool(mirparallel);
        end
        parallel = 1;
    catch
        warning('mirparallel cannot open matlabpool because it is already used.');
    end
end

if parallel
    %   The evaluation is carried out for each audio file successively
    %       (or in parallel).
    y = mirevalparallel(d,a,sr,lg,w,single,ch,export);
    isstat = 0;
else
    %   The evaluation is carried out for each audio file successively.
    y = cell(1,l);
    isstat = isfield(d,'Stat');
    for i = 1:length(order)
        f = order(i);
        if l > 1
            fprintf('\n')
            display(['*** File # ',num2str(i),'/',num2str(l),': ',a{f}]);
        end
        if mirverbose
            tic
        end
        yf = mirevalaudiofile(d,a{f},sr(f),lg(f),w(:,f),{},0,f,single,'',ch);
        if mirverbose
            toc
        end
        y{f} = yf;
        
        if (mirtemporary && length(order)>1) || not(isempty(export))
            
            if 0 % Private use.
                ff = find(a{f} == '/');
                if isempty(ff)
                    ff = 0;
                else
                    ff = ff(end);
                end
                listing = dir([a{f}(1:ff),'Fabien''s annotations/',...
                               a{f}(ff+1:end-4),'*']);
                an = importdata([a{f}(1:ff),'Fabien''s annotations/',...
                                 listing.name]);
                dt1 = mirgetdata(yf.t1);
                ok1 = find(dt1 > an*.96 & dt1 < an*1.04);
                if isempty(ok1)
                    ok1 = 0;
                end
                dt2 = mirgetdata(yf.t2);
                ok2 = find(dt2 > an*.96 & dt2 < an*1.04);
                if isempty(ok2)
                    ok2 = 0;
                end
            
                if isempty(export)
                    export = 'mirtemporary.txt';
                end
                if strncmpi(export,'Separately',10)
                    filename = a{f};
                    %filename(filename == '/') = '.';
                    %filename = ['Backup/' filename];
                    filename = [filename export(11:end)];
                    %if i == 1
                    %    mkdir('Backup');
                    %end
                    mirexport(filename,yf);
                elseif i==1
                    mirexport([export,'1'],ok1,an,yf.t1);
                    mirexport([export,'2'],ok2,an,yf.t2);
                else
                    mirexport([export,'1'],ok1,an,yf.t1,'#add');
                    mirexport([export,'2'],ok2,an,yf.t2,'#add');
                end
                
            else
                if isempty(export)
                    export = 'mirtemporary.txt';
                end
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
        end
        clear yf
    end
end

v = combineaudiofile(a,isstat,y{:});


function c = combineaudiofile(filename,isstat,varargin) % Combine output from several audio files into one single
c = varargin{1};    % The (series of) input(s) related to the first audio file
if isempty(c)
    return
end
if isstruct(c)
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
        v = cell(1,nargin-2);
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


function [l w sr lg a] = evalfolder(path,s,l,w,sr,lg,a,folders)
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
            [l w sr lg a] = evalfolder([path nf],s,l,w,sr,lg,a,1);
            cd ..
        end
    else
        [di,tpi,fpi,fi,li,bi,ni] = mirread([],nf,0,1,0);
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
                w(:,l) = min(max(interval,1),di);
            else
                w(:,l) = [1;di];
            end
            sr(l) = fi;
            lg(l) = li;
            a{l} = [path ni];
        end
    end
end