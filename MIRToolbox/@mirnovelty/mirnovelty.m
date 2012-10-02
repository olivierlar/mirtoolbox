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
%       mirnovelty(...,'Horizontal') uses the 'Horizontal' option in
%           mirsimatrix instead of 'TimeLag' used here by default.
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
    
        method.type = 'String';
        method.default = 'Foote';
        method.choice = {'Foote','Lartillot'};
    option.method = method;

            K.key = {'KernelSize','Width'};
            K.type = 'Integer';
            K.default = 64;
        option.K = K;

            gran.key = 'Granul';
            gran.type = 'Integer';
            gran.default = 1;
        option.gran = gran;

            transf.type = 'String';
            transf.default = 'TimeLag';
            transf.choice = {'Horizontal','TimeLag'};
        option.transf = transf;

            normal.key = 'Normal';
            normal.type = 'Boolean';
            normal.default = 1;
            normal.when = 'After';
        option.normal = normal;
    
specif.option = option;
specif.nochunk = 1;
varargout = mirfunction(@mirnovelty,orig,varargin,nargout,specif,@init,@main);
    

function [x type] = init(x,option)
type = 'mirscalar';
if not(isamir(x,'mirnovelty'))
    if strcmpi(option.method,'Lartillot')
        option.K = Inf;
        option.transf = 'Standard';
    end
    x = mirsimatrix(x,'Distance',option.dist,'Similarity',option.sm,...
                      'Width',max(option.K),option.transf);
end


function y = main(orig,option,postoption)
if iscell(orig)
    orig = orig{1};
end
if isa(orig,'mirnovelty')
    n = orig;
else
    s = get(orig,'Data');
    res = cell(1,length(s));
    
    if strcmpi(option.method,'Foote')
        dw = get(orig,'DiagWidth');
        Ks = option.K;
        for k = 1:length(s)
            res{k} = cell(1,length(s{k}));
            for i = 1:length(Ks)
                if isnumeric(dw)
                    dwk = dw;
                else
                    dwk = dw{k};
                end
                if Ks(i)
                    cgs = min(Ks(i),dwk);
                else
                    cgs = dwk;
                end
                cg = checkergauss(cgs,option.transf)/cgs.^option.gran;
                disp('Computing convolution, please wait...')
                for z = 1:length(s{k})
                    sz = s{k}{z};
                    szma = max(max(sz));
                    szmi = min(min(sz));
                    sz = (sz-szmi)/(szma-szmi);
                    sz = 2*sz-1;
                    sz(isnan(sz)) = 0;
                    cv = convolve2(sz,cg,'same');
                    nl = size(cv,1);
                    nc = size(cv,2);
                    if nl == 0
                        warning('WARNING IN NOVELTY: No frame decomposition. The novelty score cannot be computed.');
                        res{k}{z} = [];
                    else
                        sco = cv(floor(size(cv,1)/2),:);
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
                        res{k}{z}(i,:) = sco';
                    end
                end
            end
        end
        title = 'Novelty curve';
        
    elseif strcmpi(option.method,'Lartillot')
        for k = 1:length(s)
            res{k} = cell(1,length(s{k}));
            for z = 1:length(s{k})
                sz = s{k}{z};
                l = size(sz,1);
                pr = zeros(l);
                sel = [];
                cand = [];
                for i = 1:l
                    for j = 1:l-i
                        pr(i,j) = min(sz(i+1:i+j,i));
                    end
                    dr = find(pr(i,1:end-1)>pr(i,2:end));
                    j = 1;
                    while j <= length(cand)
                        if cand(j).current == 1
                            sel(end+1).i = cand(j).i;
                            sel(end).j = cand(j).j;
                            sel(end).sim = cand(j).sim;
                            cand(j) = [];
                        else
                            idx = find(cand(j).current-1 == dr);
                            if ~isempty(idx)
                                cand(j).current = cand(j).current-1;
                                dr(idx) = [];
                                j = j+1;
                            else
                                cand(j) = [];
                            end
                        end
                    end
                    for j = 2:length(dr)
                        cand(end+1).i = i;
                        cand(end).j = dr(j);
                        cand(end).sim = pr(i,dr(j));
                        cand(end).current = j;
                    end
                end
                res{k}{z} = sel;
                for i = 1:length(sel)
                    sel(i).i
                    sel(i).j
                    sel(i).sim
                end
            end
        end
        orig = set(orig,'Novelty',res);
        title = 'Novelty diagram';
    
    elseif strcmpi(option.method,'Lartillot.old')
        for k = 1:length(s)
            res{k} = cell(1,length(s{k}));
            for z = 1:length(s{k})
                sz = s{k}{z};
                szma = max(max(sz));
                szmi = min(min(sz));
                step = (szma-szmi)/(option.gran+1);
                thr = szmi+step:step:szma-step;
                %thr(1: floor( min(min(sz))/step )) = [];
                res{k}{z} = cell(1,length(thr));
                for h = 1:length(thr)
                    res{k}{z}{h} = [];
                    i = 1;
                    term = 0;
                    while i < size(sz,1)
                        for j = i+1:size(sz,1)
                            if find( sz(i,i+1:j) < thr(h) )
                                break
                            end
                            if j==size(sz,1)
                                term = 1;
                            end
                        end
                        if ~term
                            res{k}{z}{h}(end+1) = i+j-1;
                        end
                        i = j;
                    end
                end
            end
        end
        title = 'Novelty diagram';
        
    end
    n.method = option.method;
    n = class(n,'mirnovelty',mirdata(orig));
    n = purgedata(n);
    n = set(n,'Title',title,'Data',res,'Pos',[]);
end

if strcmp(n.method,'Foote') && ...
        not(isempty(postoption)) && postoption.normal
    s = get(orig,'Data');
    for k = 1:length(s)
        for l = 1:length(s{k})
            for i = 1:size(s{k}{l},1)
                r = s{k}{l}(i,:);
                r = (r-min(r))/(max(r)-min(r));
                s{k}{l}(i,:) = r;
            end
        end
    end
    n = set(orig,'Data',s);
end

y = {n orig};


function y = checkergauss(N,transf)
hN = ceil(N/2);
if strcmpi(transf,'TimeLag')
    y = zeros(2*N,N);
    for j = 1:N
        for i = 1:2*N+1
            g = exp(-((((i-N)-(j-hN))/hN)^2 + (((j-hN)/hN)^2))*4);
            if xor(j>hN,j-hN>i-N)
                y(i,j) = -g;
            elseif j>hN+i || j-hN<i-2*N
                y(i,j) = 0;
            else
                y(i,j) = g;
            end
        end
    end
else
    y = zeros(N);
    for i = 1:N
        for j = 1:N
            g = exp(-(((i-hN)/hN)^2 + (((j-hN)/hN)^2))*4);
            if xor(j-hN>floor((i-hN)/2),j-hN>floor((hN-i)/2))
                y(i,j) = -g;
            else
                y(i,j) = g;
            end
        end
    end
end