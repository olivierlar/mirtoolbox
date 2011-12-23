function mirsave(m,f)

if nargin == 1
    f = '.mir.mid';
end

nmat = get(m,'Data');
n = get(m,'Name');

nf = length(nmat);
for k = 1:nf
    nk = n{k};
    if nf>1 || strcmp(f(1),'.')
        nk = [nk f];
    else
        nk = f;
    end
    %writemidi(nmat{k},nk,120,60);
    nmat2midi(nmat{k},nk);
end