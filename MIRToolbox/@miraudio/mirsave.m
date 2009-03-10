function mirsave(a,f)

ext = 0;    % Specified new extension
if nargin == 1
    f = '.mir';
elseif length(f)>3 && strcmpi(f(end-3:end),'.wav')
    ext = '.wav';
    if length(f)==4
        f = '.mir';
    end
elseif length(f)>2 && strcmpi(f(end-2:end),'.au')
    ext = '.au';
    if length(f)==3
        f = '.mir';
    end
end

d = get(a,'Data');
nf = length(d);
fs = get(a,'Sampling');
nb = get(a,'NBits');
nm = get(a,'Name');
for i = 1:nf
    nbi = nb{i};
    di = d{i};
    fsi = fs{i};
    nmi = nm{i};
    
    maxd = 0;
    for j = 1:length(di)
        maxd = max(max(max(abs(di{j}),[],1),[],2),maxd);
    end

    out = [];
    for j = 1:length(di)
        di{j} = di{j}./repmat(maxd,[size(di{j},1),size(di{j},2)])*.9999;
        out = [out;reshape(di{j},[],1)];
        if length(di)>1
            out = [out;rand(100,1)]*.9;
        end
    end
    
    %Let's remove the extension from the original files
    if length(nmi)>3 && strcmpi(nmi(end-3:end),'.wav')
        nmi(end-3:end) = [];
    elseif length(nmi)>2 && strcmpi(nmi(end-2:end),'.au')
        nmi(end-2:end) = [];
    end
    
    if nf>1 || strcmp(f(1),'.')
        %Let's add the new suffix
        n = [nmi f];
    else
        n = f;
    end
    
    if not(ischar(ext)) || strcmp(ext,'.wav')
        if length(n)<4 || not(strcmpi(n(end-3:end),'.wav'))
            n = [n '.wav'];
        end
        wavwrite(out,fsi,nbi,n)
    elseif strcmp(ext,'.au')
        if length(n)<3 || not(strcmpi(n(end-2:end),'.au'))
            n = [n '.au'];
        end
        auwrite(out,fsi,nbi,'linear',n)
    end
    disp([n,' saved.']);
end