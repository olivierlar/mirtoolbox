function c = max(a,b)

d = get(a,'Data');
f = cell(1,length(d));
fp = cell(1,length(d));
n = get(a,'Name');
t = get(a,'Title');
fpa = get(a,'FramePos');

if isa(b,'mirdata')
    e = get(b,'Data');
    m = get(b,'Name');
    if not(isa(a,'miraudio'))
       t = ['max(',t,',',get(b,'Title'),')'];
    end
else
    e = {{b}};
    m = {num2str(b)};
    t = ['max(',t,',',num2str(b),')'];
end
    
for i = 1:length(d)
    f{i} = cell(1,length(d{i}));
    fp{i} = cell(1,length(d{i}));
    for j = 1:length(d{i})
        ld = size(d{i}{j},1);
        le = size(e{i}{j},1);
        %if isempty(fpa{i})
        %    ia = 1;
        %    ib = 1;
        %else
        %    [unused ia ib] = intersect(round(fpa{i}{j}(2,:)*1e4),...
        %                               round(fpb{i}{j}(2,:)*1e4));
        %end
        dj = d{i}{j};%(:,ia,:);
        ej = e{i}{j};%(:,ib,:);
        if ld > le
            ej = [ej;zeros(ld-le,size(ej,2),size(ej,3))];
        elseif ld < le
            dj = [dj;zeros(le-ld,size(dj,2),size(dj,3))];
        end
        f{i}{j} = max(dj,ej);
        if isempty(fpa{i})
            fp{i} = [];
        else
            fp{i}{j} = fpa{i}{j};%(:,ia);
        end
    end
    if isa(a,'miraudio')
        n{i} = ['max(',n{i},',',m{i},')'];
    end
end
c = set(a,'Data',f,'Name',n,'Title',t,'FramePos',fp);