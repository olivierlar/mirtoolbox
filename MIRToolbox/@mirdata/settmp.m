function x = settmp(x,tmp)

idx = get(x,'TmpIdx');
tmps = get(x,'Tmp');
tmps{idx} = tmp;
x = set(x,'Tmp',tmps);