function varargout = mirbeatspectrum(orig,varargin)
%   n = mirbeatspectrum(m) evaluates the beat spectrum from a similarity matrix.
%   [n,m] = mirbeatspectrum(m) also return the similarity matrix.
%   Optional argument: 
%       mirbeatspectrum(...,'Distance',f) specifies the name of a dissimilarity
%           distance function, from those proposed in the Statistics Toolbox
%               (help pdist).
%           default value: f = 'cosine'

% ref?

        dist.key = 'Distance';
        dist.type = 'String';
        dist.default = 'cosine';
    option.dist = dist;

        meth.type = 'String';
        meth.choice = {'Diag','Autocor'};
        meth.default = 'Autocor';
    option.meth = meth;

specif.option = option;
varargout = mirfunction(@mirbeatspectrum,orig,varargin,nargout,specif,@init,@main);
    

function [x type] = init(x,option)
if not(isamir(x,'mirscalar'))
    if isamir(x,'miraudio')
        x = mirmfcc(x,'frame',.025,'s',.01,'s','Rank',8:30);
    end
    x = mirsimatrix(x,'Distance',option.dist,'Similarity');
end
type = 'mirscalar';


function y = main(orig,option,postoption)
if iscell(orig)
    orig = orig{1};
end
fp = get(orig,'FramePos');
if not(isa(orig,'mirscalar'))
    s = get(orig,'Data');
    total = cell(1,length(s));
    for k = 1:length(s)
        maxfp = find(fp{k}{1}(2,:)>4,1);
        l = min(length(s{k}),maxfp);
        fp{k}{1}(:,maxfp+1:end) = [];
        total{k}{1} = zeros(1,l);
        if strcmpi(option.meth,'Diag')
            for i = 1:l
                total{k}{1}(i) = mean(diag(s{k},i-1));
            end
        else
            for i = 1:l
                total{k}{1}(i) = mean(mean(s{k}(:,1:l-i+1).*s{k}(:,i:l)));
            end
        end
    end
else
    total = get(orig,'Data');
end
n = mirscalar(orig,'Data',total,'FramePos',fp,'Title','Beat Spectrum'); 
y = {n orig};