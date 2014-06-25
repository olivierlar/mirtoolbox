function varargout = mirharmonicity(orig,varargin)

        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [.1 .025];
    option.frame = frame;
        
specif.option = option;
     
varargout = mirfunction(@mirharmonicity,orig,varargin,nargout,specif,@init,@main);


function [i type] = init(x,option)
if isamir(x,'miraudio')
    if option.frame.length.val
        s = mirspectrum(x,'Frame',option.frame.length.val,...
                                  option.frame.length.unit,...
                                  option.frame.hop.val,...
                                  option.frame.hop.unit,...
                                  option.frame.phase.val,...
                                  option.frame.phase.unit,...
                                  option.frame.phase.atend);
    else
        s = mirspectrum(x);
    end
else
    s = x;
end
p = mirpeaks(s,'Harmonic',20,'Contrast',.01);
i = {s,p};
type = {'mirscalar','mirscalar','mirscalar','mirspectrum'};


function h = main(x,option,postoption)
s = x{1};
p = x{2};
if iscell(p)
    p = p{1};
end
m = get(s,'Magnitude');
f = get(s,'Frequency');
pf = get(p,'TrackPos');
he = cell(1,length(m));
ie = cell(1,length(m));
hi = cell(1,length(m));
for h = 1:length(m)
    he{h} = cell(1,length(m{h}));
    ie{h} = cell(1,length(m{h}));
    hi{h} = cell(1,length(m{h}));
    for i = 1:length(m{h})
        mi = m{h}{i};
        fi = f{h}{i};
        pfi = pf{h}{i}{1};
        he{h}{i} = zeros(1,size(mi,2),size(mi,3));
        ie{h}{i} = zeros(1,size(mi,2),size(mi,3));
        hi{h}{i} = zeros(1,size(mi,2),size(mi,3));
        for j = 1:size(mi,3)
            for k = 1:size(mi,2)
                mk = mi(:,k,j);
                fk = fi(:,k,j);
                pfk = sort(pfi(:,k));
                z = zeros(2,length(pfk));
                for l = 1:length(pfk)
                    if isnan(pfk(l))
                        continue
                    end
                    f1 = find(diff(mk(pfk(l):-1:1)) >= 0,1);
                    f2 = find(diff(mk(pfk(l):end)) >= 0,1);
                    z(1,l) = pfk(l) - f1 + 1;
                    z(2,l) = pfk(l) + f2 - 1;
                end
                z(:,~z(1,:)) = [];
                hark = 0;
                inhk = sum(mk(1:z(1)-1).^2);
                for l = 1:size(z)-1
                    hark = hark + sum(mk(z(1,l):z(2,l)).^2);
                    inhk = inhk + sum(mk(z(2,l)+1:z(1,l+1)-1).^2);
                end
                inhk = inhk + sum(mk(z(end)+1:end).^2);
                he{h}{i}(1,k,j) = hark;
                ie{h}{i}(1,k,j) = inhk;
                hi{h}{i}(1,k,j) = hark / (hark + inhk);
            end
        end
    end
end
he = mirscalar(s,'Data',he,'Title','Harmonic Energy');
ie = mirscalar(s,'Data',ie,'Title','Inharmonic Energy');
hi = mirscalar(s,'Data',hi,'Title','Harmonicity');
h = {hi he ie p};