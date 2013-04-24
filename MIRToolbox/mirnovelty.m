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

        K.key = {'KernelSize','Width'};
        K.type = 'Integer';
        K.default = NaN;
    option.K = K;
    
        gran.key = 'Granul';
        gran.type = 'Integer';
        gran.default = 1;
    option.gran = gran;

        transf.type = 'String';
        transf.default = 'TimeLag';
        transf.choice = {'Horizontal','TimeLag'};
    option.transf = transf;

        flux.key = 'Flux';
        flux.type = 'Boolean';
        flux.default = 0;
    option.flux = flux;

        new.key = 'New';
        new.type = 'Boolean';
        new.default = 0;
    option.new = new;
    
        half.key = 'Half';
        half.type = 'Boolean';
        half.default = 0;
    option.half = half;

        cluster.key = 'Cluster';
        cluster.type = 'Boolean';
        cluster.default = 0;
    option.cluster = cluster;

        normal.key = 'Normal';
        normal.type = 'Boolean';
        normal.default = 1;
        normal.when = 'After';
    option.normal = normal;
    
        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [.05 1];
    option.frame = frame;

specif.option = option;
specif.nochunk = 1;
varargout = mirfunction(@mirnovelty,orig,varargin,nargout,specif,@init,@main);
    

function [x type] = init(x,option)
type = 'mirscalar';
if not(isamir(x,'mirscalar') && strcmp(get(x,'Title'),'Novelty'))
    if isnan(option.K)
        if option.cluster || option.flux || option.new
            option.K = Inf; %150;
        else
            option.K = 16; %64;
        end
    end
    if option.cluster || option.flux || option.new
        option.transf = 'Standard';
        option.half = 1;
    end
    x = mirsimatrix(x,'Distance',option.dist,'Similarity',option.sm,...
                      'Width',max(option.K),option.transf,...
                      'Half',option.half,'Cluster',option.cluster,...
                      'Frame',option.frame.length.val,option.frame.length.unit,...
                              option.frame.hop.val,option.frame.hop.unit,...
                              option.frame.phase.val,option.frame.phase.unit);
end


function y = main(orig,option,postoption)
if iscell(orig)
    orig = orig{1};
end

if option.new
    s = get(orig,'Data');
    %dw = 0; %40;
    for i = 1:length(s)
        for j = 1:length(s{i})
            %sij0 = s{i}{j}(:);
            %sij0(sij0 == 1) = [];
            %Msij = max(sij0);
            %msij = min(sij0);
            %k0 = 1;
            for k = 3:size(s{i}{j},2)
                %delt = min(k-2,dw+1);
                %dist = NaN(1,delt);
                %for l = 1:delt
                l = 1;
                    bloc = s{i}{j}(k-l-1:-1:1,1:k-l);
                    stop = find(s{i}{j}(k-l-1:-1:1,k) ...
                                > min(s{i}{j}(k-l-1:-1:1,k-l),...
                                      max(mean(bloc,2) - 2 * std(bloc,0,2),...
                                          min(bloc,[],2))));
                    if length(stop) < 2
                        stop = k-l;
                    else
                        stop(end+1) = stop(end)+1;
                        stopp = find(diff(stop(1:end-1)) == 1 & diff(stop,2) == 0,1);
                        if isempty(stopp)
                            stop = k-l;
                        else
                            stop = stop(stopp);
                        end
                    end
                    if ~isempty(stop)
                        %dist(l) = ... stop * ...
                          dist = sum(diff([s{i}{j}(k-l-stop+1:k-l-1,k),...
                                    s{i}{j}(k-l-stop+1:k-l-1,k-l)],...
                                   1,2));% ...
                                      %.*((s{i}{j}(k-l-stop:k-l-1,k-1) - msij)...
                                      %   /(Msij-msij)).^1); ...
                                  %/ max(k-1, 100);
                    end
                    %if l == 1
                    %    dist0 = dist(l);
                    %end
                    %if l > 1 && dist(l) > dist(l0) && k-l >= k0
                    %    dist(l0) = dist(l);
                    %    dist(l) = 0;
                    %else
                    %    l0 = l;
                    %    if l > 1 && k-l >= k0
                    %        dist(l0) = dist0;
                    %        k0 = k-1;
                    %    end
                    %end
                %end
                %best = 0;
                %for l = 1:delt
                    %dist(l) = dist(l) + (l-1)/(dw+1);
                %    if dist(l) > best
                %        best = dist(l);
                %        bestl = l;
                %    end
                %end
                score{i}{j}(k) = max(dist,0); %best;
                %if bestl > 1
                %    score{i}{j}(k-1:k-bestl+1) = ...
                %        min(score{i}{j}(k-1:k-bestl+1),best);
                %end
            end
        end
    end
    n = mirscalar(orig,'Data',score,'Title','Novelty'); 

elseif option.flux
    fl = mirflux(orig);
    if 1
        n = mirscalar(fl,'Title','Novelty');
    else
        score = get(fl,'Data');
        dw = get(orig,'DiagWidth');
        if dw < Inf && not(isempty(postoption)) && postoption.normal
            for k = 1:length(score)
                for l = 1:length(score{k})
                    lg = length(score{k}{l});
                    for i = 1:lg
                        score{k}{l}(i) = score{k}{l}(i) * ...
                            (1 + atan( dw/min(i,dw) - 1)) / (1 + atan(dw-1));
                    end
                end
            end
        end
        n = mirscalar(fl,'Data',score,'Title','Novelty');
    end
else
        
    fp = get(orig,'FramePos');
    if not(isa(orig,'mirscalar'))
        s = get(orig,'Data');
        cl = get(orig,'Clusters');
        if isempty(cl)
            dw = get(orig,'DiagWidth');
            Ks = option.K;
            for k = 1:length(s)
                for i = 1:length(Ks)
                    if isnumeric(dw)
                        dwk = dw;
                    else
                        dwk = dw{k};
                    end
                    if ~isnan(Ks(i)) && Ks(i)
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
                            score{k}{z} = [];
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
                            score{k}{z}(i,:) = sco';
                        end
                    end
                end
            end
        else
            score = cell(1,length(s));
            for k = 1:length(s)
                score{k} = cell(1,length(s{k}));
                for z = 1:length(s{k})
                    score{k}{z} = zeros(1,size(cl{k}{z},1));
                    for i = 1:length(cl{k}{z})
                        for j = 1:length(cl{k}{z})
                            clij = cl{k}{z}(i,j);
                            if ~isnan(clij) %&& clij > .9
                                score{k}{z}(i) = max(score{k}{z}(i), ...
                                    clij^5 * j^2);
                                score{k}{z}(i+j) = max(score{k}{z}(i+j), ...
                                    clij^5 * j^2);
                            end
                        end
                    end
                    %fp{k}{z} = [fp{k}{z}(1,2:end);fp{k}{z}(2,1:end-1)];
                end
            end
        end
    else
        score = get(orig,'Data');
    end
    if not(isempty(postoption)) && postoption.normal
        for k = 1:length(score)
            for l = 1:length(score{k})
                for i = 1:size(score{k}{l},1)
                    sco = score{k}{l}(i,:);
                    sco = (sco-min(sco))/(max(sco)-min(sco));
                    score{k}{l}(i,:) = sco;
                end
            end
        end
    end
    n = mirscalar(orig,'Data',score,'Title','Novelty','FramePos',fp); 
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