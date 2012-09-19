function mirplay(a)

ph = get(a,'Phase');
sr = get(a,'Sampling');
x = get(a,'Input');
dx = get(x,'Data');
px = get(x,'Pos');
pp = get(a,'PeakPosUnit');
for i = 1:length(pp)
    for j = 1:length(pp{i})
        for k = 1:length(pp{i}{j})
            for h = 1:length(pp{i}{j}{k})
                lag = round(pp{i}{j}{k}(h)*sr{i});
                phas = ph{i}{j}{k}(h);
                pnts = phas:lag:size(px{i}{j},1);
                dx{i}{j}(pnts,k,:) = 1;
            end
        end
    end
end
x = set(x,'Data',dx');

mirplay(x);