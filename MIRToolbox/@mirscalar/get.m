function val = get(a, propName)
% GET Get properties from the MIRdata object
% and return the value

switch propName
    case 'Mode'
        val = a.mode;
    case 'Legend'
        val = a.legend;
    case 'Parameter'
        val = a.parameter;
    case 'PeakPosUnit'
        pp = get(mirdata(a),'PeakPos');
        po = get(mirdata(a),'FramePos');
        d = get(mirdata(a),'Data');
        val = cell(1,length(pp));
        if isempty(d)
            return
        end
        for k = 1:length(pp)
            val{k} = cell(1,length(pp{k}));
            if isempty(pp{k})
                nseg = 0;
            elseif iscell(pp{k}{1})
                nseg = length(pp{k});
            else
                nseg = 1;
            end
            for i = 1:nseg
                ppi = pp{k}{i};
                if isempty(po)
                    poi = (1:size(d{k}{i},2))';
                elseif iscell(po{k})
                    if isempty(po{k})
                        poi = mean(a.framepos{k}{1})';
                    elseif ischar(po{k}{1})
                        poi = (1:length(po{k}))';
                    else
                        poi = po{k}{i};
                    end
                else
                    for j = 1:size(po,3)
                        poi(:,:,j) = po{k}(:,:,j)';
                    end
                end
                val{k}{i} = cell(size(ppi));
                for h = 1:size(ppi,3)
                    for j = 1:size(ppi,2)
%                         if size(poi,3) > 1 && size(poi,1) == 1
%                             val{k}{i}{1,j,h} = ppi{1,j,h};
%                         else
                            for x = 1:length(ppi{1,j,h})
                                val{k}{i}{1,j,h}(:,x) = ...
                                    poi(:,ppi{1,j,h}(x));
                            end
%                         end
                    end
                end
            end
        end
    otherwise
        val = get(mirdata(a),propName);
end