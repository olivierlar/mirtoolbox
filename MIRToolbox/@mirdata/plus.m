function c = plus(a,b)

d = get(a,'Data');
f = cell(1,length(d));
fp = cell(1,length(d));
n = get(a,'Name');
t = get(a,'Title');
fpa = get(a,'FramePos');

if isa(b,'mirdata')
    e = get(b,'Data');
    m = get(b,'Name');
    u = get(b,'Title');
    if not(isa(a,'miraudio'))
       t = [t,' + ',get(b,'Title')];
    end
    fpb = get(b,'FramePos');
else
    e = {{b}};
    m = {num2str(b)};
    t = [t,' + ',num2str(b)];
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
            f{i}{j} = dj + [ej;zeros(ld-le,length(ia),size(e,3))];
        elseif ld < le
            f{i}{j} = [dj;zeros(le-ld,length(ib),size(d,3))] + ej;
        else
            f{i}{j} = dj + ej;
        end
        if isempty(fpa{i})
            fp{i} = [];
        else
            fp{i}{j} = fpa{i}{j};%(:,ia);
        end
    end
    if isa(a,'miraudio')
        n{i} = [n{i} '+' m{i}];
    end
end
c = set(a,'Data',f,'Name',n,'Title',t,'FramePos',fp);