function a = subsasgn(a,index,val)
% SUBSASGN Define index assignment for mirstruct objects
switch index(1).type
case '.'
    if strcmpi(index(1).subs,'tmp')
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
end