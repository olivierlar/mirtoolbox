function val = get(a, propName)
% GET Get properties from the MIRsimatrix object
% and return the value

switch propName
    case 'DiagWidth'
        val = a.diagwidth;
    case 'Half';
        val = a.half;
    case 'Graph'
        val = a.graph;
    case 'Branch'
        val = a.branch;
    case 'Warp'
        val = a.warp;
    case 'Clusters'
        val = a.clusters;
    otherwise
        val = get(mirdata(a),propName);
end