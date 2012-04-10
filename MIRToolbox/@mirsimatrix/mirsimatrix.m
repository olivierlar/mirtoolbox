function varargout = mirsimatrix(orig,varargin)
%   m = mirsimatrix(x) computes the similarity matrix resulting from the 
%       mutual comparison between each possible frame analysis in x.
%           By default, x is the spectrum of the frame decomposition.
%           But it can be any other frame analysis.
%   Optional argument:
%       mirsimatrix(...,'Distance',f) specifies the name of a distance
%           function, from those proposed in the Statistics Toolbox
%               (help pdist).
%           default value: f = 'cosine'
%       mirsimatrix(...,'Dissimilarity') return the dissimilarity matrix,
%           which is the intermediary result before the computation of the 
%           actual similarity matrix. It shows the distance between each
%           possible frame analysis in x.
%       mirsimatrix(...,'Similarity',f) indicates the function f specifying
%           the way the distance values in the dissimilarity matrix are
%           transformed into similarity values.
%           Possible values:
%               f = 'oneminus' (default value) 
%                   corresponding to f(x) = 1-x
%               f = 'exponential'
%                   corresponding to f(x)= exp(-x)
%       mirsimatrix(...,'Width',w) or more simply dissimatrix(...,w) 
%           specifies the size of the diagonal bandwidth, in samples,
%           outside which the dissimilarity will not be computed.
%           if w = inf (default value), all the matrix will be computed.
%       mirsimatrix(...,'Horizontal') rotates the matrix 45 degrees in order to
%           make the first diagonal horizontal, and to restrict on the
%           diagonal bandwidth only.
%       mirsimatrix(...,'TimeLag') transforms the (non-rotated) matrix into
%           a time-lag matrix, making the first diagonal horizontal as well
%           (corresponding to the zero-lag line).
%
%   Foote, J. & Cooper, M. (2003). Media Segmentation using Self-Similarity
% Decomposition,. In Proc. SPIE Storage and Retrieval for Multimedia
% Databases, Vol. 5021, pp. 167-75.

%       mirsimatrix(...,'Filter',10) filter along the diagonal of the matrix,
%           using a uniform moving average filter of size f. The result is
%           represented in a time-lag matrix.


        distance.key = 'Distance';
        distance.type = 'String';
        distance.default = 'cosine';
    option.distance = distance;

        simf.key = 'Similarity';
        simf.type = 'String';
        simf.default = 'oneminus';
        simf.keydefault = 'Similarity';        
        simf.when = 'After';
    option.simf = simf;

        dissim.key = 'Dissimilarity';
        dissim.type = 'Boolean';
        dissim.default = 0;
    option.dissim = dissim;
    
        K.key = 'Width';
        K.type = 'Integer';
        K.default = Inf;
    option.K = K;
    
        filt.key = 'Filter';
        filt.type = 'Integer';
        filt.default = 0;
        filt.when = 'After';
    option.filt = filt;

        view.type = 'String';
        view.default = 'Standard';
        view.choice = {'Standard','Horizontal','TimeLag'};
        view.when = 'After';
    option.view = view;
    
        warp.key = 'Warp';
        warp.type = 'Boolean';
        warp.default = 0;
        warp.when = 'After';
    option.warp = warp;
    
        arg2.position = 2;
        arg2.default = [];
    option.arg2 = arg2;
    
specif.option = option;
specif.nochunk = 1;
varargout = mirfunction(@mirsimatrix,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if isamir(x,'miraudio')
    if isframed(x)
        x = mirspectrum(x);
    else
        x = mirspectrum(x,'Frame',0.05,1);
    end
end
type = 'mirsimatrix';


function m = main(orig,option,postoption)
if iscell(orig)
    orig = orig{1};
end
if isa(orig,'mirsimatrix')
    d = get(orig,'Data');
    for k = 1:length(d)
        if iscell(option) && isfield(option,'K') && option.K < orig.diagwidth
            nl = size(d{k},1);
            if strcmp(orig.view,'h')
                dl = floor((nl - option.K)/2);
                dk = d{k}(dl:nl-dl,:);
            else
                [spA,spd] = spdiags(d{k},-floor(option.K/2):floor(option.K/2));
                dk = full(spdiags(spA,spd,nl,size(d{k},2)));
            end
            d{k} = dk;
            orig.diagwidth = 2*floor(option.K/2)+1;
        end
    end
    m = set(orig,'Data',d);
elseif isempty(option.arg2)
    v = get(orig,'Data');
    d = cell(1,length(v));
    lK = 2*floor(option.K/2)+1;
    for k = 1:length(v)
        vk = v{k};
        if mirwaitbar
            handle = waitbar(0,'Computing dissimilarity matrix...');
        else
            handle = 0;
        end
        if not(iscell(vk))
            vk = {vk};
        end    
        for z = 1:length(vk)
            vz = vk{z};
            ll = size(vz,1);
            l = size(vz,2);
            nc = size(vz,3);
            if ll==1 && nc>1
                vz = squeeze(vz)';
                ll = nc;
                nc = 1;
            end
            nd = size(vz,4);
            if not(isempty(postoption)) && ...
                strcmpi(postoption.view,'TimeLag')
                if isinf(lK)
                    lK = l;
                end
                dk{z} = NaN(lK,l,nc);
            else
                dk{z} = NaN(l,l,nc);
            end
            for g = 1:nc
                if nd == 1
                    vv = vz;
                else
                    vv = zeros(ll*nd,l);
                    for h = 1:nd
                        if iscell(vz)
                            for i = 1:ll
                                for j = 1:l
                                    vj = vz{i,j,g,h};
                                    if isempty(vj)
                                        vv((h-1)*ll+i,j) = NaN;
                                    else
                                        vv((h-1)*ll+i,j) = vj;
                                    end
                                end
                            end
                        else
                            vv((h-1)*ll+1:h*ll,:) = vz(:,:,g,h);
                        end
                    end
                end
                if isinf(option.K) && not(strcmpi(postoption.view,'TimeLag'))
                    try
                        manually = 0;
                        dk{z}(:,:,g) = squareform(pdist(vv',option.distance));
                    catch
                        manually = 1;
                    end
                else
                    manually = 1;
                end
                if manually
                    disf = str2func(option.distance);
                    if strcmpi(option.distance,'cosine')
                        for i = 1:l
                            vv(:,i) = vv(:,i)/norm(vv(:,i));
                        end
                    end
                    if not(isempty(postoption)) && ...
                            strcmpi(postoption.view,'TimeLag')
                        hK = ceil(lK/2);
                        for i = 1:l
                            if mirwaitbar && (mod(i,100) == 1 || i == l)
                                waitbar(i/l,handle);
                            end
                            ij = min(i+lK-1,l);
                            dkij = disf(vv(:,i),vv(:,i:ij));
                            for j = 0:ij-i
                                if hK-j>0
                                    dk{z}(hK-j,i,g) = dkij(j+1);   
                                end
                                if hK+j<=lK
                                    dk{z}(hK+j,i+j,g) = dkij(j+1);
                                end
                            end
                        end
                    else
                        for i = 1:l
                            if mirwaitbar && (mod(i,100) == 1 || i == l)
                                waitbar(i/l,handle);
                            end
                            j = min(i+lK-1,l);
                            dkij = disf(vv(:,i),vv(:,i:j));
                            dk{z}(i,i:j,g) = dkij;
                            dk{z}(i:j,i,g) = dkij';
                        end
                    end
                end
            end
        end
        d{k} = dk;
        if handle
            delete(handle)
        end
    end
    m.diagwidth = lK;
    if not(isempty(postoption)) && strcmpi(postoption.view,'TimeLag')
        m.view = 'l';
    else
        m.view = 's';
    end
    m.similarity = 0;
    m.graph = {};
    m.branch = {};
    m.warp = [];
    m = class(m,'mirsimatrix',mirdata(orig));
    m = purgedata(m);
    m = set(m,'Title','Dissimilarity matrix');
    m = set(m,'Data',d,'Pos',[]);
else
    v1 = get(orig,'Data');
    v2 = get(option.arg2,'Data');
    n1 = get(orig,'Name');
    n2 = get(option.arg2,'Name');
    fp1 = get(orig,'FramePos');
    fp2 = get(option.arg2,'FramePos');
    v1 = v1{1}{1};
    v2 = v2{1}{1};
    nf1 = size(v1,2);
    nf2 = size(v2,2);
    nd = size(v1,4);
    if nd>1
        l1 = size(v1,1);
        vv = zeros(l1*nd,nf1);
        for h = 1:nd
            vv((h-1)*l1+1:h*l1,:) = v1(:,:,1,h);
        end
        v1 = vv;
        l2 = size(v2,1);
        vv = zeros(l2*nd,nf2);
        for h = 1:nd
            vv((h-1)*l2+1:h*l2,:) = v2(:,:,1,h);
        end
        v2 = vv;
        clear vv
    end
    d = NaN(nf1,nf2);
    disf = str2func(option.distance);
    if strcmpi(option.distance,'cosine')
        for i = 1:nf1
            v1(:,i) = v1(:,i)/norm(v1(:,i));
        end
        for i = 1:nf2
            v2(:,i) = v2(:,i)/norm(v2(:,i));
        end
    end
    for i = 1:nf1
        %if mirwaitbar && (mod(i,100) == 1 || i == nf1)
        %    waitbar(i/nf1,handle);
        %end
        d(i,:) = disf(v1(:,i),v2);
    end
    d = {{d}};
    m.diagwidth = NaN;
    m.view = 's';
    m.similarity = 0;
    m.graph = {};
    m.branch = {};
    m.warp = [];
    m = class(m,'mirsimatrix',mirdata(orig));
    m = purgedata(m);
    m = set(m,'Title','Dissimilarity matrix','Data',d,'Pos',[],...
              'Name',{n1{1},n2{1}},'FramePos',{fp1{1},fp2{1}});
end
lK = option.K;
if not(isempty(postoption))
    if strcmpi(m.view,'s') && isempty(option.arg2)
        if strcmpi(postoption.view,'Horizontal')
            for k = 1:length(d)
                for z = 1:length(d{k})
                    d{k}{z} = rotatesim(d{k}{z},m.diagwidth);
                    if lK < m.diagwidth
                        W = size(d{k}{z},1);
                        hW = ceil(W/2);
                        hK = floor(lK/2);
                        d{k}{z} = d{k}{z}(hW-hK:hW+hK,:);
                        m.diagwidth = lK;
                    end
                end
            end
            m = set(m,'Data',d);
            m.view = 'h';
        elseif strcmpi(postoption.view,'TimeLag') || postoption.filt
            W = m.diagwidth;
            for k = 1:length(d)
                for z = 1:length(d{k})
                    if isinf(W)
                        dz = NaN(size(d{k}{z}));
                    else
                        dz = NaN(W,size(d{k}{z},2));
                    end
                    for l = 1:size(dz,1)
                        dz(l,l:end) = diag(d{k}{z},l-1)';
                    end
                    if lK < m.diagwidth
                        W = size(dz,1);
                        hK = ceil(lK/2);
                        if size(dz,1)>hK
                            dz = dz(1:hK,:);
                        end
                        m.diagwidth = lK;
                    end
                    d{k}{z}= dz;
                end
            end
            m = set(m,'Data',d);
            m.view = 'l';
        end
    end
    if ischar(postoption.simf)
        if strcmpi(postoption.simf,'Similarity')
            if not(isequal(m.similarity,NaN))
                option.dissim = 0;
            end
            postoption.simf = 'oneminus';
        end
        if isequal(m.similarity,0) && isstruct(option) ...
                && isfield(option,'dissim') && not(option.dissim)
            simf = str2func(postoption.simf);
            for k = 1:length(d)
                for z = 1:length(d{k})
                     d{k}{z} = simf(d{k}{z});
                end
            end
            m.similarity = postoption.simf;
            m = set(m,'Title','Similarity matrix','Data',d);
        elseif length(m.similarity) == 1 && isnan(m.similarity) ...
                && option.dissim
            m.similarity = 0;
        end
    end
    if postoption.filt
        fp = get(m,'FramePos');
        for k = 1:length(d)
            for z = 1:length(d{k})
                dz = filter(ones(postoption.filt,1),1,d{k}{z});
                d{k}{z} = dz(postoption.filt:end,1:end-postoption.filt+1);
                fp{k}{z} = [fp{k}{z}(1,1:end-postoption.filt+1);...
                            fp{k}{z}(1,postoption.filt:end)];
            end
        end
        m = set(m,'Data',d,'FramePos',fp);
    end
    if postoption.warp
        dz = 1 - (d{1}{1} - min(min(d{1}{1}))) / (max(max(d{1}{1})) - min(min(d{1}{1})));
        dz(dz>.33) = Inf;
        bests = zeros(size(dz)+1);
        bests(1,2:end) = Inf;
        bests(2:end,1) = Inf;
        paths = cell(size(dz)+1);
        for i = 1:size(dz,1)
            for j = 1:size(dz,2)
                [best,index] = min([bests(i,j),...
                                    bests(i+1,j),...
                                    bests(i,j+1)]);
                bests(i+1,j+1) = best + dz(i,j);
                if ~isinf(bests(i+1,j+1))
                    switch index
                        case 1
                            path = paths{i,j};
                        case 2
                            path = paths{i+1,j};
                        case 3
                            path = paths{i,j+1};
                    end
                    paths{i+1,j+1} = [path;i j];
                end
            end
        end
        scod = repmat(((1:size(dz,1))/size(dz,1))',[1 size(dz,2)]) .* ...
               repmat(((1:size(dz,2))/size(dz,2)),[size(dz,1) 1]);
        [unused sord] = sort(scod(:),'descend');
        paths(1,:) = [];
        paths(:,1) = [];
        for i = 1:length(sord)
            if ~isempty(paths{sord(i)})
                m = set(m,'Warp',paths{sord(i)});
                break
            end
        end
    end
end


function S = rotatesim(d,K)
if length(d) == 1;
    S = d;
else
    K = min(K,size(d,1)*2+1);
    lK = floor(K/2);
    S = NaN(K,size(d,2),size(d,3));
    for k = 1:size(d,3)
        for j = -lK:lK
            S(lK+j+1,:,k) = [NaN(1,floor(abs(j)/2)) diag(d(:,:,k),j)' ...
                                                    NaN(1,ceil(abs(j)/2))];
        end
    end
end

function d = cosine(r,s)
d = 1-r'*s;
%nr = sqrt(r'*r);
%ns = sqrt(s'*s);
%if or(nr == 0, ns == 0);
%    d = 1;
%else
%    d = 1 - r'*s/nr/ns;
%end


function d = KL(x,y)
% Kullback-Leibler distance
if size(x,4)>1
    x(end+1:2*end,:,:,1) = x(:,:,:,2);
    x(:,:,:,2) = [];
end
if size(y,4)>1
    y(end+1:2*end,:,:,1) = y(:,:,:,2);
    y(:,:,:,2) = [];
end
m1 = mean(x);
m2 = mean(y);
S1 = cov(x);
S2 = cov(y);
d = (trace(S1/S2)+trace(S2/S1)+(m1-m2)'*inv(S1+S2)*(m1-m2))/2 - size(S1,1);
    

function s = exponential(d)
    s = exp(-d);
    
    
function s = oneminus(d)
    s = 1-d;