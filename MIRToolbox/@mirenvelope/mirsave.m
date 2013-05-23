function mirsave(e,f)

ext = 0;    % Specified new extension
if nargin == 1
    f = '.envelope.mir';
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

d = get(e,'Data');
nf = length(d);
fs = get(e,'Sampling');
nm = get(e,'Name');
pp = get(e,'PeakPosUnit');
for i = 1:nf
    di = d{i}{1};
    nmi = nm{i};
    di = resample(di,11025,round(fs{i}));
    di = rand(size(di)).*di;
    di = di/max(max(max(abs(di))))*.9999;
    di = reshape(di,[],1);
    
    if ~isempty(pp) && ~isempty(pp{i}{1})
        pi = pp{i}{1}{1};
        d2i = zeros(length(di),1);
        for h = 1:length(pi)
            d2i(round(pi(h)*11025)) = 1;
        end
        di = di/10 + d2i;
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
        wavwrite(di,11025,32,n)
    elseif strcmp(ext,'.au')
        if length(n)<3 || not(strcmpi(n(end-2:end),'.au'))
            n = [n '.au'];
        end
        auwrite(di,11025,32,'linear',n)
    end
    disp([n,' saved.']);
end