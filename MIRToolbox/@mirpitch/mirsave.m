function mirsave(a,f)

d = get(a,'Data');
amp = get(a,'Amplitude');
mef = get(a,'Mean');
n = get(a,'Name');
fp = get(a,'FramePos');

if nargin == 1
    f = '.mirpitch.wav';
end

nf = length(d);
for k = 1:nf
    dk = d{k};
    if not(iscell(dk))
        dk = {dk};
    end
    nk = n{k};
    out = [];
    for l = 1:size(dk{1},3)
        for i = 1:length(dk)
            di = dk{i};
            synth = zeros(ceil((fp{k}{i}(end)-fp{k}{i}(1))*44100)+1,1);
            if isempty(mef)
                transcribed = 0;
            else
                transcribed = ~isempty(mef{k}{i});
            end
            if transcribed
                stp = get(a,'Start');
                enp = get(a,'End');
                mef = mef{k}{i}{1};
                stp = stp{k}{i}{1};
                enp = enp{k}{i}{1};
                for j = 1:length(mef)
                    fj = 2^(mef(j)/1200);
                    k1 = floor((fp{k}{i}(1,stp(j))-fp{k}{i}(1))*44100)+1;
                    k2 = floor((fp{k}{i}(2,enp(j))-fp{k}{i}(1))*44100)+1;
                    synth(k1:k2) = synth(k1:k2) ...
                        + sum(sin(2*pi*fj*(0:k2-k1)'/44100),2) ...
                                .*hann(k2-k1+1);
                end
            else
                ampi = amp{k}{i};
                for j = 1:size(di,2)
                    if iscell(di)
                        dj = di{j};
                    else
                        dj = di(:,j);
                    end
                    dj(isnan(dj)) = 0;
                    ampj = zeros(size(dj));
                    if iscell(ampi)
                        ampj(1:size(ampi{j})) = ampi{j};
                    else
                        ampj(1:size(ampi(:,j))) = ampi(:,j);
                    end
                    if not(isempty(dj))
                        k1 = floor((fp{k}{i}(1,j)-fp{k}{i}(1))*44100)+1;
                        k2 = floor((fp{k}{i}(2,j)-fp{k}{i}(1))*44100)+1;
                        ampj = repmat(ampj,1,k2-k1+1);
                        synth(k1:k2) = synth(k1:k2) ...
                            + sum(ampj.*sin(2*pi*dj*(0:k2-k1)/44100),1)' ...
                                    .*hann(k2-k1+1);
                    end
                end
            end
            out = [out;synth];
            if size(dk{1},3)>1
                out = [out;rand(1,10)];
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