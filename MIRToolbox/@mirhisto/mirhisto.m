function varargout = mirhisto(x,varargin)
%   h = mirhisto(x) constructs the histogram from x. The elements of x are
%       binned into equally spaced containers.
%   Optional argument:
%       mirhisto(...,'Number',n): specifies the number of containers.
%           Default value : n = 10.
%       mirhisto(...,'Ampli'): adds the amplitude of the elements,instead of
%           simply counting then.


        n.key = 'Number';
        n.type = 'Integer';
        n.default = 10;
    option.n = n;
    
        a.key = 'Ampli';
        a.type = 'Boolean';
        a.default = 0;
    option.a = a;
    
specif.option = option;


varargout = mirfunction(@mirhisto,x,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
type = 'mirhisto';


function h = main(x,option,postoption)
if iscell(x)
    x = x{1};
end
d = get(x,'Data');
fp = get(x,'FramePos');
%disp('Computing histogram...')
ddd = cell(1,length(d));
bbb = cell(1,length(d));
for i = 1:length(d)
    di = d{i}{1};  % To be generalized for segmented data
    if iscell(di)
        mx = -Inf;
        mn = Inf;
        nc = size(di,2);
        for k = 1:nc
            dk = di{k};
            if size(dk,4) == 2
                dk(end+1:end*2,:,:,1) = dk(:,:,:,2);
                dk(:,:,:,2) = [];
            end
            mxk = max(dk);
            mnk = min(dk);
            if mxk > mx
                mx = mxk;
            end
            if mnk < mn
                mn = mnk;
            end
        end
        if isinf(mx) || isinf(mx)
            b = [];
            dd = [];
        else
            dd = zeros(1,option.n);
            if mn == mx
                b(1,:) = mn-ceil(option.n/2) : mn+floor(option.n/2);
            else
                b(1,:) = mn : (mx-mn)/option.n : mx;
            end
            for k = 1:nc
                dk = di{k};
                for j = 1:option.n
                    found = find(and(dk>=b(1,j),dk<=b(1,j+1)));
                    if option.a
                        dd(1,j) = dd(1,j) + sum(dk(found));
                    else
                        dd(1,j) = dd(1,j) + length(found);
                    end
                end
            end
        end
    else
        if isa(x,'mirscalar')
            di = permute(di,[3 2 1]);
        end
        if size(di,4) == 2
            di(end+1:end*2,:,:,1) = di(:,:,:,2);
            di(:,:,:,2) = [];
        end
        nl = size(di,1);
        nc = size(di,2);
        np = size(di,3);
        dd = zeros(1,option.n,np);
        for l = 1:np
            mx = max(max(di(:,:,l),[],1),[],2);
            mn = min(min(di(:,:,l),[],1),[],2);
            b(l,:) = mn:(mx-mn)/option.n:mx;
            for k = 1:nc
                dk = di(:,k,l);
                for j = 1:option.n
                    found = (find(and(dk>=b(l,j),dk<=b(l,j+1))));
                    if option.a
                        dd(1,j,l) = dd(1,j,l) + sum(dk(found));
                    else
                        dd(1,j,l) = dd(1,j,l) + length(found);
                    end
                end
            end
        end
    end
    ddd{i} = ipermute(dd,[3 2 1]);
    bbb{i}(:,:,1) = b(:,1:end-1);
    bbb{i}(:,:,2) = b(:,2:end);
    fp{i} = {fp{i}{1}([1 end])'};
end
x = set(x,'FramePos',fp);
h = class(struct,'mirhisto',mirdata(x));
h = purgedata(h);
h = set(h,'Bins',bbb,'Weight',ddd);