function a = subsasgn(a,index,val)
% SUBSASGN Define index assignment for mirstruct objects
switch index(1).type
case '.'
    if strcmpi(index(1).subs,'tmp')
        fields = a.fields;
        data = a.data;
        
        if isa(val,'mirdata') || (iscell(val) && isa(val{1},'mirdata'))
            % If the 'tmp' data turns to be an actual evaluated data,
            % the mirstruct object is transformed into a simple struct.
            a = struct;
            for i = 1:length(fields)
                a.(fields{i}) = data{i};
            end
            a.tmp.(index(2).subs) = val;
            return
        end
        
        if isa(val,'mirdesign')
            val = set(val,'Stored',{index.subs});
        end
        if length(index)>2
            if strcmpi(index(3).type,'{}')
                isubs = index(3).subs;
                if length(isubs)>1
                    a.tmp.(index(2).subs){isubs{1},isubs{2}} = val;
                else
                    a.tmp.(index(2).subs){isubs{1}} = val;
                end
            end
        else
            a.tmp.(index(2).subs) = val;
        end
        aa = struct;
        aa.fields = fields;
        aa.data = data;
        aa.tmp = a.tmp;
        aa.stat = a.stat;
        a = class(aa,'mirstruct',val);
        return
    end
    [is,id] = ismember(index(1).subs,a.fields);
    if not(is)
        a.fields{end+1} = index(1).subs;
        a.data{end+1} = [];
        id = length(a.fields);
    end
    if length(index) == 1
        a.data{id} = val;
    else
        a.data{id} = subsasgn(a.data{id},index(2:end),val);
    end
    if get(val,'NoChunk') && isframed(a)
        a = set(a,'FrameDontChunk',1);
        % Frame-decomposed flowchart where one dependent variable requires
        % a complete computation. Should not therefore be evaluated
        % chunk after chunk.
    end
end