function varargout = mirnovelty(orig,varargin)
%   n = mirnovelty(m) evaluates the novelty score from a similarity matrix.
%   [n,m] = mirnovelty(m) also return the similarity matrix.
%   Optional argument: 
%       mirnovelty(...,'Distance',f) specifies the name of a dissimilarity
%           distance function, from those proposed in the Statistics Toolbox
%               (help pdist).
%           default value: f = 'cosine'
%       mirnovelty(...,'Similarity',f) specifies the name of a similarity
%           measure function. This function is applied to the result of the
%           distance function. cf. mirsimatrix 
%           default value: f = 'exponential'
%               corresponding to f(x) = exp(-x)
%       mirnovelty(...,'KernelSize',s) or more simply mirnovelty(...,s) 
%           specifies the length of the gaussian kernel, in samples.
%           default value: s = 64.
%       mirnovelty(...,'Normal',0) does not normalize the novelty curve
%           between the values 0 and 1.
%
%   Foote, J. & Cooper, M. (2003). Media Segmentation using Self-Similarity
% Decomposition,. In Proc. SPIE Storage and Retrieval for Multimedia
% Databases, Vol. 5021, pp. 167-75.

        dist.key = 'Distance';
        dist.type = 'String';
        dist.default = 'cosine';
    option.dist = dist;

        sm.key = {'Measure','Similarity'};
        sm.type = 'String';
        sm.default = 'exponential';
    option.sm = sm;

        K.key = {'KernelSize','Width'};
        K.type = 'Integer';
        K.default = 64;
    option.K = K;

        normal.key = 'Normal';
        normal.type = 'Boolean';
        normal.default = 1;
        normal.when = 'After';
    option.normal = normal;
    
specif.option = option;
specif.combineframes = @combineframes;
specif.nochunk = 1;
varargout = mirfunction(@mirnovelty,orig,varargin,nargout,specif,@init,@main);
    

function [x type] = init(x,option)
type = 'mirscalar';
if not(isamir(x,'mirscalar') && strcmp(get(x,'Title'),'Novelty'))
    x = mirsimatrix(x,'Distance',option.dist,'Similarity',option.sm,...
                      'Width',option.K,'Horizontal');
end
if isa(x,'mirdesign')
    x = set(x,'Overlap',ceil(option.K));
end


function y = main(orig,option,postoption)
if iscell(orig)
    orig = orig{1};
end
if not(isa(orig,'mirscalar'))
    s = get(orig,'Data');
    dw = get(orig,'DiagWidth');
    for k = 1:length(s)
        if isnumeric(dw)
            dwk = dw;
        else
            dwk = dw{k};
        end
        if option.K
            cgs = min(option.K,dwk);
        else
            cgs = dwk;
        end
        hgs = floor(cgs/2);
        cg = checkergauss(cgs);
        disp('Computing convolution, please wait...')
        for z = 1:length(s{k})
            sz = s{k}{z};
            szm = max(max(sz));
            for i = find(isnan(sz))
                sz(i) = szm;
            end
            cv = convolve2(sz,cg,'same');
            nl = size(cv,1);
            nc = size(cv,2);
            if nl == 0
                warning('WARNING IN NOVELTY: No frame decomposition. The novelty score cannot be computed.');
                score{k}{z} = [];
            else
                sco = cv(ceil(nl/2),:);
                incr = find(diff(sco)>=0);
                if not(isempty(incr))
                    decr = find(diff(sco)<=0);
                    sco(1:incr(1)-1) = NaN(1,incr(1)-1);
                    if not(isempty(decr))
                        sco(decr(end)+1:end) = NaN(1,length(sco)-decr(end));
                    end
                    incr = find(diff(sco)>=0);
                    sco2 = sco;
                    if not(isempty(incr))
                        sco2 = sco2(1:incr(end)+1);
                    end
                    decr = find(diff(sco)<=0);
                    if not(isempty(decr)) && decr(1)>2
                        sco2 = sco2(decr(1)-1:end);
                    end
                    mins = min(sco2);
                    rang = find(sco>= mins);
                    if not(isempty(rang))
                        sco(1:rang(1)-1) = NaN(1,rang(1)-1);
                        sco(rang(end)+1:end) = NaN(1,length(sco)-rang(end));
                    end
                end
                score{k}{z} = sco;
            end
        end
    end
else
    score = get(orig,'Data');
end
if not(isempty(postoption)) && postoption.normal
    for k = 1:length(score)
        for l = 1:length(score{k})
            sco = score{k}{l};
            sco = (sco-min(sco))/(max(sco)-min(sco));
            score{k}{l} = sco;
        end
    end
end
n = mirscalar(orig,'Data',score,'Title','Novelty'); 
y = {n orig};


function old = combineframes(old,new)
if not(iscell(old))
    old = {old};
end
if not(iscell(new))
    new = {new};
end
for var = 1:length(new)
    ov = old{var};
    nv = new{var};
    ofp = get(ov,'FramePos');
    ofp = ofp{1}{1};
    nfp = get(nv,'FramePos');
    nfp = nfp{1}{1};
    od = get(ov,'Data');
    od = od{1}{1};
    onan = find(isnan(od));
    od(onan) = [];
    ofp(:,onan) = [];
    nd = get(nv,'Data');
    nd = nd{1}{1};
    nnan = find(isnan(nd));
    nd(nnan) = [];
    nfp(:,nnan) = [];
    [unused omatch nmatch] = intersect(ofp(1,:),nfp(1,:));
    if isempty(omatch)
        ov = set(ov,'FramePos',{{[ofp nfp]}},'Data',{{[od nd]}});
    else
        lm = length(omatch);
        ov = set(ov,'FramePos',{{[ofp(:,1:omatch(1)-1) nfp]}},...
            'Data',{{[od(1:omatch(1)-1),...
                      (od(omatch).*(lm:-1:1) + nd(nmatch).*(1:lm))/(lm+1),...
                      nd(nmatch(end)+1:end)]}});
    end
    old{var} = ov;
end

function y = checkergauss(N)
hN = ceil(N/2);
y = zeros(N);
for j = 1:N
    for i = 1:N
        g = exp(-(((i-hN)/hN)^2 + (((j-hN)/hN)^2))*4);
        if j>hN && xor(j>i,j>N-i)
            y(i,j) = -g;
        else
            y(i,j) = g;
        end
    end
end
y(hN,hN) = 1;