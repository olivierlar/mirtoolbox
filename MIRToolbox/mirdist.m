function d = mirdist(x,y,dist)
%   d = mirdist(x,y) evaluates the distance between x and y.
%   Optional argument:
%       d = mirdist(x,y,f) specifies distance function.
%           Default value: 'Cosine'

clx = get(x,'Clusters');
if isempty(clx{1})
    px = get(x,'PeakPos');
    if not(iscell(px)) || isempty(px{1}) || ...
            not(iscell(px{1})) || not(iscell(px{1}{1}))
        if nargin < 3
            dist = 'Cosine';
        end

        d = get(x,'Data');
        d = d{1}{1};
        if iscell(d)
            d = d{1};
        end

        e = get(y,'Data');
        dt = cell(1,length(e));
        for h = 1:length(e)
            ee = e{h}{1};
            if iscell(ee)
                ee = ee{1};
            end
            if length(d)<length(ee)
                ee = ee(1:length(d));
                dd = d;
            else
                dd = d(1:length(ee));
            end
            if length(dd) == 1
                dt{h}{1} = abs(dd-ee);
            else
                dt{h}{1} = pdist([dd(:)';ee(:)'],dist);
            end
        end
    else
        sig = pi/4;
        dx = get(x,'Data');
        nx = length(px{1}{1}{1});
        cx = mean(px{1}{1}{1}/length(dx{1}{1}));
        c
        dy = get(y,'Data');
        py = get(y,'PeakPos');
        dt = cell(1,length(py));
        for h = 1:length(py)
            ny = length(py{h}{1}{1});
            cy = mean(py{h}{1}{1}/length(dy{h}{1}));
            dt{h}{1} = sqrt((nx*cos(sig*cx)-ny*cos(sig*cy))^2 ...
                           +(nx*sin(sig*cx)-ny*sin(sig*cy))^2);
        end
    end
else
    cly = get(y,'Clusters');
    dt = cell(1,length(cly));
    for h = 1:length(cly)
        cost = zeros(length(clx{1}.weight),length(cly{h}.weight));
        for i = 1:length(clx{1}.weight)
            for j = 1:length(cly{h}.weight)
                covx = clx{1}.covar(:,i);
                covy = cly{h}.covar(:,j);
                mux = clx{1}.centr(:,i);
                muy = cly{h}.centr(:,j);
                cost(i,j) = sum(covx./covy + covy./covx + ...
                        (mux-muy).^2.*(1./covx + 1./covy) - 2);
            end
        end
        dt{h}{1} = emd_wrapper(cost,clx{1}.weight,cly{h}.weight);
    end
end
d = mirscalar(y,'Data',dt,'Title',[get(y,'Title'),' Distance'],...
                'Name',get(x,'Name'),'Name2',get(y,'Name'));