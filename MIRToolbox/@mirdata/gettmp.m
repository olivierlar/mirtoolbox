function [tmp y] = gettmp(x,y)

idx = get(x,'TmpIdx')+1;
tmps = get(x,'Tmp');
if idx > length(tmps)
    tmp = [];
else
    tmp = tmps{idx};
end
if nargin<2
    y = x;
end
y = set(y,'Tmp',tmps,'TmpIdx',idx);