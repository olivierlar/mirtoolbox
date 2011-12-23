function mirsave(a,f)

d = get(a,'Data');
n = get(a,'Name');
fp = get(a,'FramePos');

if nargin == 1
    t = get(a,'Title');
    f = ['.' t '.mir'];
end

nf = length(d);
for k = 1:nf
    dk = d{k};
    nk = n{k};
    if not(iscell(dk))
        dk = {dk};
    end
    out = [];
    for l = 1:size(dk{1},3)
        for i = 1:length(dk)
            di = dk{i};
            synth = zeros(1,ceil((fp{k}{i}(end)-fp{k}{i}(1))*44100)+1);
            for j = 1:size(di,2)
                if iscell(di)
                    dj = di{j};
                else
                    dj = di(j);
                end
                if not(isempty(dj))
                    k1 = floor((fp{k}{i}(1,j)-fp{k}{i}(1))*44100)+1;
                    k2 = floor((fp{k}{i}(2,j)-fp{k}{i}(1))*44100)+1;
                    synth(k1:k2) = synth(k1:k2) ...
                        + sum(sin(2*pi*dj*(0:k2-k1)/44100),1).*hann(k2-k1+1)';
                end
            end
            out = [out synth];
            if size(dk{1},3)>1
                out = [out rand(1,10)];
            end
        end
    end
    fout = miraudio(out,44100);
    
    if nf>1 || strcmp(f(1),'.')
        nk = [nk f];
    else
        nk = f;
    end
    mirsave(fout,nk);
end