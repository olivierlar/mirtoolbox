function v = mirevalaudiofile(d,file,sampling,lg,size,struc,istmp,index,single,name,ch)
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
        res = mirevalaudiofile(field,file,sampling,lg,size,struc,istmp,index,...
                                                     single,fieldname,ch);
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
        v{j} = mirevalaudiofile(d{j},file,sampling,lg,size,struc,istmp,index,...
                                       single,[name,num2str(j)],ch);
    end
elseif isa(d,'mirstruct') && isempty(get(d,'Argin'))
    mirerror('MIRSTRUCT','You should always use tmp fields when using mirstruct. Else, just use struct.');
elseif get(d,'SeparateChannels')
    v = cell(1,ch);
    for i = 1:ch
        d = set(d,'File',file,'Sampling',sampling,'Length',lg,'Size',size,...
                  'Eval',1,'Index',index,'Struct',struc,'Channel',i);
        % For that particular file or this particular feature, let's begin the
        % actual evaluation process.
        v{i} = evaleach(d,single,name);    
        % evaleach performs a top-down traversal of the design flowchart.
    end
    v = combinechannels(v);
else
    d = set(d,'File',file,'Sampling',sampling,'Length',lg,'Size',size,...
              'Eval',1,'Index',index,'Struct',struc);
    dl = get(d,'FrameLength');
    if length(dl)>1
        v = cell(1,length(dl));
        for i = 1:length(dl)
            d = set(d,'Scale',i);
            v{i} = evaleach(d,single,name);
        end
        v = combinescales(v);
    else
        % For that particular file or this particular feature, let's begin the
        % actual evaluation process.
        v = evaleach(d,single,name);    
        % evaleach performs a top-down traversal of the design flowchart.
    end
end


function y = combinechannels(c)
y = c{1};
v = get(y,'Data');
for h = 2:length(c)
    d = get(c{h},'Data');
    for i = 1:length(d)
        if isa(y,'mirmidi')
            d{i}(:,3) = h;
            v{i} = sortrows([v{i};d{i}]);
        else
            for j = 1:length(d{i})
                v{i}{j}(:,:,h) = d{i}{j};
            end
        end
    end
end
y = set(y,'Data',v);


function y = combinescales(s)
y = s{1};
fp = get(y{1},'FramePos');
fp = fp{1};
for j = 1:length(y)
    v = get(y{j},'Data');
    for h = 2:length(s)
        d = get(s{h}{j},'Data');
        for i = 1:length(d)
            v{i}{h} = d{i}{1};
        end
        if j == 1
            fph = get(s{h}{j},'FramePos');
            fp{h} = fph{1}{1};
        end
    end
    y{j} = set(y{j},'Data',v,'FramePos',{fp});
end