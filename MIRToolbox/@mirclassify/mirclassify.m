function c = mirclassify(x,varargin)
%   Optional argument:
%       mirclassify(...,'Nearest') uses the minimum distance strategy.
%           (by default)
%       mirclassify(...,'Nearest',k) uses the k-nearest-neighbour strategy.
%           Default value: k = 1, corresponding to the minimum distance
%               strategy.
%       mirclassify(...,'GMM',ng) uses a gaussian mixture model. Each class is
%           modeled by at most ng gaussians.
%           Default value: ng = 1.
%           Additionnally, the type of mixture model can be specified,
%           using the set of value proposed in the gmm function: i.e.,
%           'spherical','diag','full' (default value) and 'ppca'.
%               (cf. help gmm)
%           Requires the Netlab toolbox.

if not(iscell(x))
    x = {x};
end
n = get(x{1},'Name');

l = length(n);    % Number of training samples
lab = cell(1,l);
allabs = struct;
for i = 1:l
    sl = strfind(n{i},'/');
    if isempty(sl)
        mirerror('mirclassify','There should not be audio files in the main folder. They should be all in subfolders corresponding to the different classes.');
    end
    lab{i} = n{i}(1:sl(end)-1);
    if i == 1
        allabs.name = lab{1};
        allabs.idx = 1;
        continue
    end
    [test li] = ismember(lab{i},{allabs.name});
    if test
        allabs(li).idx(end+1) = i;
    else
        allabs(end+1).name = lab{i};
        allabs(end).idx = i;
    end
end

[k,ncentres,covartype,kmiter,emiter,d,mahl] = scanargin(varargin);
nfolds = length(allabs(1).idx);
v = [];                        % Preprocessed training vectors
mn = cell(1,length(x));
sd = cell(1,length(x));
for i = 1:length(x)
    if isnumeric(x{i})
        d = cell(1,size(x{i},2));
        for j = 1:size(x{i},2)
            d{j} = x{i}(:,j);
        end
    else
        d = get(x{i},'Data');
    end
    
    [v mn{i} sd{i}] = integrate(v,d,l);
    
    if 0 %isa(t{i},'scalar')
        m = mode(x{i});
        if not(isempty(m))
            v = integrate(v,m,l);
        end
    end
end
mahl = cov(v');

for i = 1:length(allabs)
    s = RandStream('mt19937ar','Seed','shuffle');
    p = randperm(s,length(allabs(i).idx));
    allabs(i).idx = allabs(i).idx(p);
end

for h = 1:nfolds
    va = [];
    vt = [];
    idxh = zeros(1,length(allabs));
    rlab = cell(1,length(allabs));
    rlabt = cell(1,length(allabs)*(nfolds-1));
    for i = 1:length(allabs)
        idxi = allabs(i).idx;
        va = [va v(:,allabs(i).idx(h))];
        for j = (nfolds-1)*(i-1)+(1:nfolds-1)
            rlabt{j} = allabs(i).name;
        end
        idxi(h) = [];
        vt = [vt v(:,idxi)];  
        rlab{i} = allabs(i).name;
    end

    if k                % k-Nearest Neighbour
        c.nbparam = size(vt,2);
        for l = 1:size(va,2)
            [sv,idx] = sort(distance(va(:,l),vt,mahl));
            labs = cell(0); % Class labels
            founds = [];    % Number of found elements in each class
            for i = idx(1:k)
                labi = rlabt{i};
                found = 0;
                for j = 1:length(labs)
                    if isequal(labi,labs{j})
                        found = j;
                    end
                end
                if found
                    founds(found) = founds(found)+1;
                else
                    labs{end+1} = labi;
                    founds(end+1) = 1;
                end
            end
            [b ib] = max(founds);
            c.classes{h,l} = labs{ib};
        end
    elseif ncentres     % Gaussian Mixture Model
        labs = cell(0);    % Class labels
        founds = cell(0);  % Elements associated to each label.
        for i = 1:size(vt,2)
            labi = rlabt{i};
            found = 0;
            for j = 1:length(labs)
                if isequal(labi,labs{j})
                    founds{j}(end+1) = i;
                    found = 1;
                end
            end
            if not(found)
                labs{end+1} = labi;
                founds{end+1} = i;
            end
        end
        options      = zeros(1, 18);
        options(2:3) = 1e-4;
        options(4)   = 1e-6;
        options(16)  = 1e-8;
        options(17)  = 0.1;
        options(1)   = 0; %Prints out error values, -1 else
        c.nbparam = 0;
        OK = 0;
        while not(OK)
            OK = 1;
            for i = 1:length(labs)
                options(14)  = kmiter;
                try
                    mix{i} = gmm(size(vt,1),ncentres,covartype);
                catch
                    error('ERROR IN CLASSIFY: Netlab toolbox not installed.');
                end
                mix{i} = netlabgmminit(mix{i},vt(:,founds{i})',options);
                options(5)   = 1;
                options(14)  = emiter;
                try
                    mix{i} = gmmem(mix{i},vt(:,founds{i})',options);
                    c.nbparam = c.nbparam + ...
                        length(mix{i}.centres(:)) + length(mix{i}.covars(:));
                catch
                    %err = lasterr;
                    %warning('WARNING IN CLASSIFY: Problem when calling GMMEM:');
                    %disp(err);
                    %disp('Let us try again...');
                    OK = 0;
                end
            end    
        end
        pr = zeros(size(va,2),length(labs));
        for i = 1:length(labs)
            prior = length(founds{i})/size(vt,2);
            pr(:,i) = prior * gmmprob(mix{i},va');
            %c.post{i} = gmmpost(mix{i},va');
        end
        [mm ib] = max(pr');
        for i = 1:size(va,2)
            c.classes{h,i} = labs{ib(i)};
        end
    end
    if isempty(rlab)
        c.correct = NaN;
    else
        correct = 0;
        for i = 1:length(rlab)
            if isequal(c.classes{h,i},rlab{i})
                correct = correct + 1;
            end
        end
        c.correct(h) = correct / length(rlab);
    end   
end

c = class(c,'mirclassify');

disp('Expected classes:')
rlab

disp('Classification results:')
c.classes

if isnan(c.correct)
    disp('No label has been associated to the test set. Correct classification rate cannot be computed.');
else
    disp(['Correct classification rate: ',num2str(mean(c.correct),2)]);
end

figure
if k
    v0 = zeros(size(v,1),0);
    for i = 1:length(allabs)
        v0 = [v0 v(:,allabs(i).idx)];
    end
    ltick = nfolds;
elseif ncentres
    v0 = zeros(size(v,1),length(mix)*ncentres);
    for i = 1:length(mix)
        for j = 1:ncentres
            v0(:,ncentres*(i-1)+j) = mix{i}.centres(j,:)';
        end
    end
    ltick = ncentres;
end
imagesc(v0)
set(gca,'XTick',ltick+.5:ltick:size(v,2)+.5)
set(gca,'XTickLabel','')
set(gca,'YTick',[])
grid on
for i = 1:length(allabs)
    text(ltick*(i-1)+1,1,allabs(i).name);
end

function [vt m s] = integrate(vt,v,lvt,m,s)
% lvt is the number of samples
vtl = [];
for l = 1:lvt
    vl = v{l};
    if iscell(vl)
        vl = vl{1};
    end
    if iscell(vl)
        vl = vl{1};
    end
    if size(vl,2) > 1
        mirerror('MIRCLASSIFY','The analytic features guiding the classification should not be frame-decomposed.');
    end
    vtl(:,l) = vl;
end

if nargin<4
    m = mean(vtl,2);
    s = std(vtl,0,2);
end

dnom = repmat(s,[1 size(vtl,2)]);
dnom = dnom + (dnom == 0);  % In order to avoid division by 0
vtl = (vtl - repmat(m,[1 size(vtl,2)])) ./ dnom;

vt(end+1:end+size(vtl,1),:) = vtl;


function [k,ncentres,covartype,kmiter,emiter,d,mahl] = scanargin(v)
k = 1;
d = 0;
i = 1;
ncentres = 0;
covartype = 'full';
kmiter = 10;
emiter = 100;
mahl = 1;
while i <= length(v)
    arg = v{i};
    if ischar(arg) && strcmpi(arg,'Nearest')
        k = 1;
        if length(v)>i && isnumeric(v{i+1})
            i = i+1;
            k = v{i};
        end
    elseif ischar(arg) && strcmpi(arg,'GMM')
        k = 0;
        ncentres = 1;
        if length(v)>i
            if isnumeric(v{i+1})
                i = i+1;
                ncentres = v{i};
                if length(v)>i && ischar(v{i+1})
                    i = i+1;
                    covartype = v{i};
                end
            elseif ischar(v{i+1})
                i = i+1;
                covartype = v{i};
                if length(v)>i && isnumeric(v{i+1})
                    i = i+1;
                    ncentres = v{i};
                end
            end                
        end
    elseif isnumeric(arg)
        k = v{i};
    else
        error('ERROR IN MIRCLASSIFY: Syntax error. See help mirclassify.');
    end    
    i = i+1;
end


function y = distance(a,t,mahl)

for i = 1:size(t,2)
    if det(mahl) > 0  % more generally, uses cond
        lham = inv(mahl);
    else
        lham = pinv(mahl);
    end
    y(i) = sqrt((a - t(:,i))'*lham*(a - t(:,i)));        
end
%y = sqrt(sum(repmat(a,[1,size(t,2)])-t,1).^2);