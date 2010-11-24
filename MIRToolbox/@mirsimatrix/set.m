function mm = set(m,varargin)
% SET Set properties for the MIRsimatrix object
% and return the updated object

propertyArgIn = varargin;
graph = m.graph;
branch = m.branch;
d = mirdata(m);
while length(propertyArgIn) >= 2,
   prop = propertyArgIn{1};
   val = propertyArgIn{2};
   propertyArgIn = propertyArgIn(3:end);
   switch prop
       case 'Graph'
           graph = val;
       case 'Branch'
           branch = val;
       otherwise
           d = set(d,prop,val);
   end
end
mm.diagwidth = m.diagwidth;
mm.view = m.view;
mm.similarity = m.similarity;
mm.graph = graph;
mm.branch = branch;
mm = class(mm,'mirsimatrix',d);