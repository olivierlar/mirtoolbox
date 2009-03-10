function b = subsref(a,index)
% Define field name indexing for mirstruct objects
switch index(1).type
case '.'
    if strcmpi(index(1).subs,'tmp')
        if length(index)== 1
            b = [];
        else
            if length(index)>2
                if strcmpi(index(3).type,'{}')
                    isubs = index(3).subs;
                    if length(isubs)>1
                        b = a.tmp.(index(2).subs){isubs{1},isubs{2}};
                    else
                        b = a.tmp.(index(2).subs){isubs{1}};
                    end
                end
            else
                b = a.tmp.(index(2).subs);
            end
        end
        return
    end
    [is,id] = ismember(index(1).subs,a.fields);
    if length(index) == 1
        b = a.data{id};
    else
        b = subsref(a.data{id},index(2:end));
    end
end