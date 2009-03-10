function mirplay(q,varargin)

rank = 1:length(q.name);
if nargin>2 && strcmpi(varargin{1},'Sequence')
    rank = intersect(rank,varargin{2});
end
disp(['Query: ',q.query.name{1}]);
mirplay(q.query.name{1});
disp('********')
disp('Retrieval:')
for i = rank
    disp([num2str(i),'. ',q.name{i},' (dist = ',num2str(q.dist(i)),')']);
    mirplay(q.name{i});
end