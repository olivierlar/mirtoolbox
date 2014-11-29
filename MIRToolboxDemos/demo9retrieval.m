function demo9retrieval(query)

if nargin<1
    query = 'vivaldi.wav';
end

%%
% 1. Timbre

cc0 = mirmfcc(query,'frame',.025,'s',.01,'s','Rank',8:30)
cl0 = mircluster(cc0,16)

cc = mirmfcc('Folder','frame',.025,'s',.01,'s','Rank',8:30);
cl = mircluster(cc,16);

d1 = mirdist(cl0,cl)
mirquery(cl0,cl)
mirplay(ans,'Sequence',1:3)

pause

%%
% 2. Rhythm

[bs0 sm0] = mirbeatspectrum(query)
bs = mirbeatspectrum('Folder');

d2 = mirdist(bs0,bs)
mirquery(bs0,bs)
mirplay(ans,'Sequence',1:3)

    % Variant
    
    [tp0,ac0] = mirtempo(query) 
    [tp,ac] = mirtempo('Folder');
    ac0 = purgedata(ac0);
    ac = purgedata(ac);
    d2bis = mirdist(ac0,ac)
    mirquery(ac0,ac)
    mirplay(ans,'Sequence',1:3)

pause
    
%%
% 3. Structure

pk0 = mirpeaks(mirnovelty(query))
pk = mirpeaks(mirnovelty('Folder'));

d3 = mirdist(pk0,pk)
mirquery(pk0,pk)
mirplay(ans,'Sequence',1:3)

pause

%%
% 4. Combination

d = d1*.6 + d2*.3 + d3*.1
mirquery(d)
mirplay(ans,'Sequence',1:3)

dbis = d1*.6 + d2bis*.3+ d3*.1
mirquery(dbis)
mirplay(ans,'Sequence',1:3)