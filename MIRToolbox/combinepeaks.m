function c = combinepeaks(p,v,thr)
% dedicated function for (Klapuri, 99) that creates a curve made of burst
% at position of peaks p and with amplitude related to peaks p2.

d = get(v,'Data');
pp = get(p,'PeakPos');
pv = get(v,'PeakPos');
sr = get(v,'Sampling');
l = length(d);
empty = cell(1,l);
for i = 1:l
    thr = round(thr*sr{i});
    empty{i} = cell(1,length(d{i}));
    for h = 1:length(d{i})
        dih = zeros(size(d{i}{h}));
        for l = 1:size(pp{i}{h},3)
            for k = 1:size(pp{i}{h},2)
                j = 1;
                ppkl = pp{i}{h}{1,k,l};
                pvkl = pv{i}{h}{1,k,l};
                while j < length(ppkl)
                    if ppkl(j+1)-ppkl(j) < thr
                        decreas = d{i}{h}(pvkl(j+1),k,l) ...
                                  < d{i}{h}(pvkl(j),k,l);
                        ppkl(j+decreas) = [];
                        pvkl(j+decreas) = [];
                    else
                        j = j+1;
                    end
                end
                dih(ppkl,k,l) = d{i}{h}(pvkl,k,l);
            end
        end
        d{i}{h} = dih;
    end
end
c = set(v,'Data',d,'PeakPos',empty,'PeakVal',empty,...
          'PeakPrecisePos',{},'PeakPreciseVal',{},'PeakMode',empty,...
          'Tmp',[]);