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
%           diagonal bandwith only.
%
%       mirsimatrix(M,r) creates a mirsimatrix similarity matrix based on
%           the Matlab square matrix M, of frame rate r (in Hz.)
%               By default r = 20 Hz.
%           mirsimatrix(M,r,'Dissimilarity') creates instead a mirsimatrix
%               dissimilarity matrix.
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
    
        rate.type = 'Integer';
        rate.position = 2;
        rate.default = 20;
    option.rate = rate;
    
specif.option = option;
specif.nochunk = 1;
varargout = mirfunction(@mirsimatrix,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if isnumeric(x)
    m.diagwidth = Inf;
    m.view = 's';
    m.similarity = NaN;
    m = class(m,'mirsimatrix',mirdata);
    m = set(m,'Title','Dissimilarity matrix');
    fp = repmat(((1:size(x,1))-.5)/option.rate,[2,1]);
    x = set(m,'Data',{x},'Pos',[],...
              'FramePos',{{fp}},'Name',{inputname(1)});
else
    if not(isamir(x,'mirsimatrix'))
        if (isamir(x,'miraudio'))
            if isframed(x)
                x = mirspectrum(x);
            else
                x = mirspectrum(x,'Frame',0.05,1);
            end
        end
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
                dl = (nl - option.K)/2;
                dk = d{k}(ceil(dl):nl-floor(dl),:);
            else
                [spA,spd] = spdiags(d{k},-ceil(option.K/2):ceil(option.K/2));
                dk = full(spdiags(spA,spd,nl,size(d{k},2)));
            end
            d{k} = dk;
            orig.diagwidth = option.K;
        end
    end
    m = set(orig,'Data',d);
else
    v = get(orig,'Data');
    d = cell(1,length(v));
    for k = 1:length(v)
        vk = v{k};
        hK = floor(option.K/2);
        if mirwaitbar
            handle = waitbar(0,'Computing dissimilarity matrix...');
        else
            handle = 0;
        end
        if 0 %iscell(vk)
            try
                vk = cell2mat(vk);
            end
        end
        if 0 %iscell(vk) %&& length(vk) > 1 %%% ATTENTION KL!!<<<<<<<<<<<<
            l = length(vk);
            dk = NaN(l,l);
            for i = 1:l
                if handle
                    waitbar(i/l,handle);
                end
                for j = max(1,i-hK):min(l,i+hK)
                    dk(i,j) = KL(vk{i},vk{j});
                end
            end
        else
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
                    strcmpi(postoption.view,'Horizontal')
                    dk{z} = NaN(option.K,l,nc);
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
                                tic
                                vv((h-1)*ll+1:h*ll,:) = vz(:,:,g,h);
                                toc
                            end
                        end
                    end
                    if isinf(option.K)  
                        try
                            manually = 0;
                            dk{z}(:,:,g) = squareform(pdist(vv',option.distance));
                            if option.K < Inf
                                [spA,spd] = spdiags(dk{z},...
                                    -ceil(option.K/2):ceil(option.K/2));
                                dk{z}(:,:,g) = full(spdiags(spA,spd,size(dk,1),size(dk{z},2)));
                            end
                        catch
                            %err = lasterror;
                            %warning(err.message)
                            %disp('Statistics Toolbox does not seem to be
                            %installed. Recompute the distance matrix manually.');
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
                                strcmpi(postoption.view,'Horizontal')
                            for i = 1:l
                                if mirwaitbar && (mod(i,100) == 1 || i == l)
                                    waitbar(i/l,handle);
                                end
                                lK = floor(option.K/2);
                                ij = min(i+lK-1,l);
                                dkij = disf(vv(:,i),vv(:,i:ij));
                                for j = 1:ij-i+1
                                    dk{z}(lK+j-1,i+j-1,g) = dkij(j);
                                    dk{z}(lK-j+1,i+j-1,g) = dkij(j);
                                end
                            end
                        else
                            for i = 1:l
                                if mirwaitbar && (mod(i,100) == 1 || i == l)
                                    waitbar(i/l,handle);
                                end
                                j = min(i+option.K-1,l);
                                dkij = disf(vv(:,i),vv(:,i:j));
                                dk{z}(i,i:j,g) = dkij;
                                dk{z}(i:j,i,g) = dkij';
                            end
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
    m.diagwidth = option.K;
    if not(isinf(option.K)) && not(isempty(postoption)) && ...
            strcmpi(postoption.view,'Horizontal')
        m.view = 'h';
    else
        m.view = 's';
    end
    m.similarity = 0;
    m = class(m,'mirsimatrix',mirdata(orig));
    m = purgedata(m);
    m = set(m,'Title','Dissimilarity matrix');
    m = set(m,'Data',d,'Pos',[]);
end
if not(isempty(postoption))
    if strcmpi(m.view,'s')
        if strcmpi(postoption.view,'Horizontal')
            for k = 1:length(d)
                for z = 1:length(d{k})
                    d{k}{z} = rotatesim(d{k}{z},m.diagwidth);
                end
            end
            m = set(m,'Data',d);
            m.view = 'h';
        elseif strcmpi(postoption.view,'TimeLag') || postoption.filt
            for k = 1:length(d)
                for z = 1:length(d{k})
                    dz = NaN(size(d{k}{z}));
                    for l = 1:size(dk,1)
                        dz(1:end-l+1,l) = diag(d{k}{z},l-1);
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
end


function S = rotatesim(d,K)
if length(d) == 1;
    S = d;
else
    K = min(K,size(d,1)*2);
    lK = floor(K/2);
    S = NaN(K,size(d,2),size(d,3));
    for k = 1:size(d,3)
        for j = -lK:lK
            S(lK+1+j,:,k) = [NaN(1,floor(abs(j)/2)) diag(d(:,:,k),j)' ...
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