function c = mirclassify(a,da,t,dt,varargin)
%   c = mirclassify(test,features_test,train,features_train) classifies the
%       audio sequence(s) contained in the audio object test, along the
%       analytic feature(s) features_test, following the supervised
%       learning of a training set defined by the audio object train and
%       the corresponding analytic feature(s) features_train.
%           * The analytic feature(s) features_test should *not* be frame 
%               decomposed. Frame-decomposed data should first be 
%               summarized, using for instance mirmean or mirstd.
%           * Multiple analytic features have to be grouped into one array 
%               of cells.
%       You can also integrate your own arrays of numbers computed outside
%           MIRtoolbox as part of the features. These arrays should be 
%           given as matrices where each successive column is the analysis 
%           of each successive file.
%   Example:
%       mirclassify(test, mfcc(test), train, mfcc(train))
%       mirclassify(test, {mfcc(test), centroid(test)}, ...
%                train, {mfcc(train), centroid(train)})
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

lab = get(t,'Label');
c.labtraining = lab;
rlab = get(a,'Label');
c.labtest = rlab;
[k,ncentres,covartype,kmiter,emiter,d,mahl] = scanargin(varargin);
disp('Classifying...')

%% Training phase
if not(iscell(dt))
    dt = {dt};
end
lvt = length(get(t,'Data'));    % Number of training samples
vt = [];                        % Preprocessed training vectors
mn = cell(1,length(dt));
sd = cell(1,length(dt));
for i = 1:length(dt)
    if isnumeric(dt{i})
        d = cell(1,size(dt{i},2));
        for j = 1:size(dt{i},2)
            d{j} = dt{i}(:,j);
        end
    else
        d = get(dt{i},'Data');
    end
    
    [vt mn{i} sd{i}] = integrate(vt,d,lvt);
    
    if 0 %isa(dt{i},'scalar')
        m = mode(dt{i});
        if not(isempty(m))
            vt = integrate(vt,m,lvt);
        end
    end
end
c.training = vt;
dim = size(vt,1);

%% Test phase
if not(iscell(da))
    da = {da};
end
lva = length(get(a,'Data'));    % Number of test samples
va = [];                        % Preprocessed test vectors
for i = 1:length(da)
    if isnumeric(da{i})
        d = cell(1,size(da{i},2));
        for j = 1:size(da{i},2)
            d{j} = da{i}(:,j);
        end
    else
        d = get(da{i},'Data');
    end
    
    va = integrate(va,d,lva,mn{i},sd{i});
    
    if 0 %isa(da{i},'scalar')
        m = mode(da{i});
        if not(isempty(m))
            va = integrate(va,m,lva,m{i},s{i});
        end
    end
end
c.test = va;

c.nbobs = lvt;
totva = [vt va];
mahl = cov(totva');
if k                % k-Nearest Neighbour
    c.nbparam = lvt;
    for l = 1:lva
        [sv,idx] = sort(distance(va(:,l),vt,d,mahl));
        labs = cell(0); % Class labels
        founds = [];    % Number of found elements in each class
        for i = idx(1:k)
            labi = lab{i};
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
        c.classes{l} = labs{ib};
    end
elseif ncentres     % Gaussian Mixture Model
    labs = cell(0);    % Class labels
    founds = cell(0);  % Elements associated to each label.
    for i = 1:lvt
        labi = lab{i};
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
                mix{i} = gmm(dim,ncentres,covartype);
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
                err = lasterr;
                warning('WARNING IN CLASSIFY: Problem when calling GMMEM:');
                disp(err);
                disp('Let us try again...');
                OK = 0;
            end
        end    
    end
    pr = zeros(lva,length(labs));
    for i = 1:length(labs)
        prior = length(founds{i})/lvt;
        pr(:,i) = prior * gmmprob(mix{i},va');
        %c.post{i} = gmmpost(mix{i},va');
    end
    [mm ib] = max(pr');
    for i = 1:lva
        c.classes{i} = labs{ib(i)};
    end
end
if isempty(rlab)
    c.correct = NaN;
else
    correct = 0;
    for i = 1:lva
        if isequal(c.classes{i},rlab{i})
            correct = correct + 1;
        end
    end
    c.correct = correct / lva;
end        
c = class(c,'mirclassify');


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


function y = distance(a,t,d,mahl)

for i = 1:size(t,2)
    if det(mahl) > 0  % more generally, uses cond
        lham = inv(mahl);
    else
        lham = pinv(mahl);
    end
    y(i) = sqrt((a - t(:,i))'*lham*(a - t(:,i)));        
end
%y = sqrt(sum(repmat(a,[1,size(t,2)])-t,1).^2);