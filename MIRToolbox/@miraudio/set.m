function aa = set(a,varargin)
% SET Set properties for the miraudio object
% and return the updated object

t = mirtemporal(a);
t = set(t,'Title',get(a,'Title'),'Abs',get(a,'Abs'),'Ord',get(a,'Ord'),...
        varargin{:});
aa.fresh = a.fresh;
aa = class(aa,'miraudio',t);