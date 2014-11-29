function [a d] = mircluster(a,varargin)
%   c = mircluster(a,f) clusters the segments in the audio sequence(s) 
%       contained in the audio object a, along the analytic feature(s) 
%       f, using the k-means strategy. Multiple analytic features have to
%       be grouped into one array of cells.
%       Example:
%           sg = mirsegment(a);
%           mircluster(sg, mirmfcc(sg))
%           mircluster(sg, {mirmfcc(sg), mircentroid(sg)})
%   c = mircluster(d) clusters the frame-decomposed data d into groups
%       using K-means clustering.
%       Example:
%           cc = mirmfcc(a,'Frame');
%           mircluster(cc)
%   Optional argument:
%       mircluster(...,n) indicates the maximal number of clusters.
%           Default value: n = 2.
%       mircluster(...,'Runs',r) indicates the maximal number of runs.
%           Default value: r = 5.
%
%   Requires SOM Toolbox (included in the MIRtoolbox distribution).


        nruns.key = 'Runs';
        nruns.type = 'Integer';
        nruns.default = 5;
    option.nruns = nruns;
     
        nclust.position = 2;
        nclust.type = 'Integer';
        nclust.default = 2;
    option.nclust = nclust;

specif.option = option;

specif.nochunk = 1;
    
d = a;
if isa(a,'mirdesign')
    if not(get(a,'Eval'))
        % During bottom-up construction of the general design

        [unused option] = miroptions(@mircluster,a,specif,varargin);
        type = get(a,'Type');
        a = mirdesign(@mircluster,a,option,{},struct,type);  
        a = set(a,'NoChunk',1);
    else
        % During top-down evaluation initiation
        
        e = evaleach(a);
        if iscell(e)
            e = e{1};
        end
        a = mircluster(e,varargin{:});
    end
else
    if not(isa(a,'mirdata'))
        mirerror('mircluster','The input should be either frame- or segment-decomposed.');
    end
   
    if isempty(varargin) || (not(isa(varargin{1},'mirdata') || ...
                             (iscell(varargin{1}) && ...
                              isa(varargin{1}{1},'mirdata'))))
        % mircluster version for frame-decomposed data:
        % frames are clustered into groups using K-means clustering.
        [unused option] = miroptions(@mircluster,a,specif,varargin);
        da = get(a,'Data');
        lva = length(da); % Number of audio files in the audio object.
        c = cell(1,lva);
        display('Clustering frames...');
        if mirwaitbar
            handle = waitbar(0,'Clustering frames...');
        else
            handle = 0;
        end
        for j = 1:lva           % For each audio file,...
            va = [];            % Data transmitted to the kmeans_cluster function.
            v = da{j};
            if iscell(v)
                v = uncell(v,-Inf); %v{1};
            end
            if size(v,4)>1
                v(end+1:2*end,:,:,1) = v(:,:,:,2);
                v(:,:,:,2) = [];
            end
            % Standardization
            %stv = std(v,0,2);
            %stv(find(stv == 0)) = 1;
            va(end+1:end+size(v,1),:,:) = v;%...
                %(v - repmat(mean(v,2),[1 size(v,2) ])) ...
                %./ repmat(stv,[1 size(v,2) ]);
            if isa(a,'mirscalar')
                m = get(a,'Mode');
                if not(isempty(m))
                    m = m{j};
                    val = [];
                    for l = 1:nseg
                        vl = m{l};
                        if iscell(vl)
                            vl = vl{1};
                        end
                        val(:,l) = vl;
                    end
                    stv = std(val,0,2);
                    stv(find(stv == 0)) = 1;
                    va(end+1:end+size(val,1),:) = ...
                        (val - repmat(mean(val,2),[1 size(val,2) ])) ...
                            ./ repmat(stv,[1 size(val,2) ]);
                end
            end
            if size(va,3)>1
                mel = 1;
                va = reshape(va,size(va,2),size(va,3))';
            else
                mel = 0;
            end
            [cc, p, err, ind] = kmeans_clusters(va',option.nclust,option.nruns);
            [minind select] = min(ind);
            c{j}.centr = cc{select}';
            c{j}.index = p{select};
            c{j}.weight = zeros(1,size(cc{select},1));
            c{j}.covar = zeros(size(cc{select}'));
            ii = 1;
            for i = 1:size(c{j}.centr,2)
                clus = va(:,c{j}.index == ii);
                if isempty(clus)
                    higher = find(c{j}.index > ii);
                    c{j}.index(higher) = c{j}.index(higher)-1;
                    c{j}.centr(:,ii) = [];
                    c{j}.weight(ii) = [];
                    c{j}.covar(:,ii) = [];
                else
                    c{j}.weight(ii) = size(clus,2)/size(va,2);
                    c{j}.covar(:,ii) = mean((clus'-ones(1,size(clus,1))*c{j}.centr(:,ii)).^2);
                    ii = ii+1;
                end
            end
            if isa(handle,'matlab.ui.Figure')
                waitbar(j/lva,handle);
            end
        end
        if isa(handle,'matlab.ui.Figure')
           delete(handle)
        end
        a = set(a,'Clusters',c);
    else
        % mircluster version for segmented audio:
        % segments are clustered into groups using K-means clustering.
        da = varargin{1};
        varargin(1) = [];
        [unused option] = miroptions(@mircluster,a,specif,varargin);
        display('Clustering segments...');
        if isa(da,'mirdata') || (iscell(da) && isa(da{1},'mirdata'))
            if not(iscell(da))
                da = {da};
            end
            vala = get(a,'Data');    % Data contained in the audio object a.
            lva = length(vala); % Number of audio files in the audio object.
            clus = cell(1,lva);
            for j = 1:lva           % For each audio file,...
                va = [];            % Data transmitted to the kmeans_cluster function.
                nseg = length(vala{j}); % Number of segments in the audio file.
                for i = 1:length(da)    % For each analytic feature,...
                    v = get(da{i},'Data');
                    v = v{j};
                    if iscell(v)
                        v = uncell(v,-Inf); %v{1};
                    end
                    val = [];           
                    if size(v,4)>1
                        v(end+1:2*end,:,:,1) = v(:,:,:,2);
                        v(:,:,:,2) = [];
                    end

                    % Standardization
                    stv = std(v,0,2);
                    stv(find(stv == 0)) = 1;
                    va(end+1:end+size(v,1),:) = ...
                        (v - repmat(mean(v,2),[1 size(v,2) ])) ...
                        ./ repmat(stv,[1 size(v,2) ]);
                    if isa(da{i},'mirscalar')
                        m = get(da{i},'Mode');
                        if not(isempty(m))
                            m = m{j};
                            val = [];
                            for l = 1:nseg
                                vl = m{l};
                                if iscell(vl)
                                    vl = vl{1};
                                end
                                val(:,l) = vl;
                            end
                            stv = std(val,0,2);
                            stv(find(stv == 0)) = 1;
                            va(end+1:end+size(val,1),:) = ...
                                (val - repmat(mean(val,2),[1 size(val,2) ])) ...
                                    ./ repmat(stv,[1 size(val,2) ]);
                        end
                    end

                end
                [cc, p, err, ind] = kmeans_clusters(va',min(option.nclust,nseg),option.nruns);
                clus{j} = p{end};
            end
            a = set(a,'Clusters',clus);
            t = get(a,'Time'); 
            fp = get(a,'FramePos'); 
            for j = 1:lva           % For each audio file,...
                aj = vala{j};
                tj = t{j};
                fpj = fp{j};
                clj = clus{j};
                k = 2;
                while k <= length(aj)
                    if clj(k) == clj(k-1)
                        aj{k-1} = [aj{k-1};aj{k}];
                        aj(k) = [];
                        tj{k-1} = [tj{k-1};tj{k}];
                        tj(k) = [];
                        fpj{k-1} = [fpj{k-1}(1);fpj{k}(2)];
                        fpj(k) = [];
                        clj(k) = [];
                        k = k-1;
                    end
                    k = k+1;
                end
                vala{j} = aj;
                t{j} = tj;
                fp{j} = fpj;
                cl{j} = clj;
            end
            a = set(a,'Data',vala,'Time',t,'FramePos',fp,'Clusters',cl);
        end
    end
end