function mirplay(e,varargin)

% mirplay method for mirsimatrix objects.

specif.option = struct;

specif.eachchunk = 'Normal';

varargout = mirfunction(@mirplay,e,varargin,nargout,specif,@init,@main);
if nargout == 0
    varargout = {};
end


function [x type] = init(x,option)
type = '';


function noargout = main(m,option,postoption)
d = get(m,'Data');
n = get(m,'Name');
fp = get(m,'FramePos');
w = get(m,'Warp');

figure
fp1 = cell2mat(fp{1});
fp2 = cell2mat(fp{2});
fp1m = (fp1(1,:,:,:)+fp1(2,:,:,:))/2;
fp2m = (fp2(1,:,:,:)+fp2(2,:,:,:))/2;
imagesc(fp2m,fp1m,d{1}{1});
hold on

paths = w{1};
bests = w{2};
bestsindex = w{3};
[i j] = find(bests);
mirverbose(0)
for k = 1:length(i)
    best = bests(i(k),j(k));
    bestindex = bestsindex(i(k),j(k));
    path = paths{real(best)+1,imag(best)+1}{bestindex};
    if path(2,end) - path(2,1) > 10 && path(1,end) - path(1,1) > 10
        if 1
            for h = 1:size(path,2)
                plot(fp2m(path(2,h)),fp1m(path(1,h)),'k+')
            end
            drawnow
            pause
            for h = 1:size(path,2)
                plot(fp2m(path(2,h)),fp1m(path(1,h)),'w+')
            end
            %a = miraudio(n{1},'Extract',fp1(1,path(1,1)),...
            %                            fp1(2,path(1,end)));
            %b = miraudio(n{1},'Extract',fp2(1,path(2,1)),...
            %                            fp2(2,path(2,end)));
            %mirplay(a);
            %mirplay(b);
        else
            d = 1;
            h = 1;
            while h < size(path,2)
                chge = find(path(d,h+1:end) > path(d,h),1);
                if isempty(chge)
                    hh = h:size(path,2);
                else
                    hh = h:h+chge-1;
                end
                plot(fp2m(path(2,hh)),fp1m(path(1,hh)),'k+');
                drawnow
                a = miraudio(n{1},'Extract',fp1(1,path(1,hh(1))),...
                                            fp1(2,path(1,hh(end))));
                b = miraudio(n{1},'Extract',fp2(1,path(2,hh(1))),...
                                            fp2(2,path(2,hh(end))));
                mirplay(a);
                mirplay(b);
                h = hh(end)+1;
                d = 3-d;
            end
        end            
    end
end
noargout = {};