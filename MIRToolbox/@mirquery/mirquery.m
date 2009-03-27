function res = mirquery(varargin)
%   r = mirquery(q,b), where
%           q is the analysis of one audio file and 
%           b is the analysis of a folder of audio files,
%               according to the same mirtoolbox feature,
%       returns the name of the audio files in the database b in an
%           increasing distance to q with respect to the chosen feature.
%   r = mirquery(d), where
%           d is the distance between one audio file and a folder of audio
%           file, according to a mirtoolbox feature,
%       returns the name of the audio files in an increasing distance d.
%
%   Optional argument:
%       mirquery(...,'Best',n) returns the name of the n closest audio
%           files.
%       mirquery(..,'Distance',d) specifies the distance to use.
%           Default value: d = 'Cosine'   (cf. mirdist)

[distfunc,nbout] = scanargin(varargin);

if nargin<2 || not(isa(varargin{2},'mirdata'))
    returnval=0;
    dist = varargin{1};
    name =  get(dist,'Name2');
    res.query.val = [];
    res.query.name = get(dist,'Name');
else
    returnval=1;
    query = varargin{1};
    base = varargin{2};
    name = get(base,'Name');
    res.query.val = mirgetdata(query);
    res.query.name = get(query,'Name');
    database = mirgetdata(base);
    dist = mirdist(query,base,distfunc);
end
datadist = mirgetdata(dist);

[ordist order] = sort(datadist);
order(isnan(ordist)) = [];
nbout = min(nbout,length(order));
res.dist = ordist(1:nbout);
if returnval
    res.val = database(order);
else
    res.val = [];
end
res.name = name(order);

res = class(res,'mirquery');



function [distfunc,nbout] = scanargin(v)
distfunc = 'Cosine';
nbout=Inf;
i = 1;
while i <= length(v)
    arg = v{i};
    if ischar(arg) && strcmpi(arg,'Distance')
        if length(v)>i && ischar(v{i+1})
            i = i+1;
            distfunc = v{i};
        end
    elseif ischar(arg) && strcmpi(arg,'Best')
        if length(v)>i && isnumeric(v{i+1})
            i = i+1;
            nbout = v{i};
        end
    %else
    %    error('ERROR IN MIRQUERY: Syntax error. See help mirquery.');
    end    
    i = i+1;
end