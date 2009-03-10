function display(q)

% MIRQUERY/DISPLAY display of mirquery

disp(['Query: ',q.query.name{1}]);
disp('********')
disp('Retrieval:')
for i = 1:length(q.dist)
    disp([num2str(i),'. ',q.name{i},' (dist = ',num2str(q.dist(i)),')']);
end