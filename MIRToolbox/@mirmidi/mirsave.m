function mirsave(m,f)

ext = 0;    % Specified new extension
if nargin == 1
    f = '.mir';
elseif length(f)>3 && strcmpi(f(end-3:end),'.mid')
    ext = '.mid';
    if length(f)==4
        f = '.mir';
    end
elseif length(f)>2 && strcmpi(f(end-2:end),'.ly')
    ext = '.ly';
    if length(f)==3
        f = '.mir';
    end
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

    if not(ischar(ext)) || strcmp(ext,'.mid')
        if length(n)<4 || not(strcmpi(n(end-3:end),'.mid'))
            n = [n '.mid'];
        end
        %writemidi(nmat{k},nk,120,60);
        nmat2midi(nmat{k},nk);
    elseif strcmp(ext,'.ly')
        if length(n)<3 || not(strcmpi(n(end-2:end),'.ly '))
            n = [n '.ly'];
        end
        lywrite(nmat{k},nk);
    end
    disp([nk,' saved.']);
end


function lywrite(nmat,filename)
fid = fopen(filename,'wt');
v = ver('MIRtoolbox');
fprintf(fid,['% LilyPond score automatically generated using MIRtoolbox version ' v.Version]);
fprintf(fid,'\\new Score \\with {\\override TimeSignature #''transparent = ##t} \n');
fprintf(fid,'\\relative c'' { \n');
fprintf(fid,'\\cadenzaOn \n');
for i = 1:size(nmat,1)
    switch mod(nmat(i,4)+2,12)
        case 0
            p = 'c';
        case 1
            p = 'cis';
        case 2
            p = 'd';
        case 3
            p = 'dis';
        case 4
            p = 'e';
        case 5
            p = 'f';
        case 6
            p = 'fis';
        case 7
            p = 'g';
        case 8
            p = 'gis';
        case 9
            p = 'a';
        case 10
            p = 'ais';
        case 11
            p = 'b';
    end
    if ~mod(i,15)
        fprintf(fid,'\\bar ""\n');
    end
    if i>1
        do = round((nmat(i,4)-nmat(i-1,4))/12);
        if do>0
            for j = 1:do
                p = [p ''''];
            end
        elseif do<0
            for j = 1:-do
                p = [p ','];
            end
        end
    end
    du = nmat(i,2);
    if du < .2
        fprintf(fid,[p '16 ']);
    else
        fprintf(fid,[p '8 ']);
    end
end
fprintf(fid,'\n } \n \n \\version "2.14.2"');