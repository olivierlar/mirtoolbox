function varargout = mirpattern(orig,varargin)
%   p = mirpattern(a)

        period.key = 'Period';
        period.type = 'Boolean';
        period.when = 'After';
        period.default = 0;
    option.period = period;
        
specif.option = option;
     
varargout = mirfunction(@mirpattern,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if not(isamir(x,'mirpattern'))
    x = mirsimatrix(x);
end
type = 'mirpattern';


function p = main(orig,option,postoption)
if not(isamir(orig,'mirpattern'))
    b = get(orig,'Branch');
    fp = get(orig,'FramePos');
    pp = get(orig,'PeakPos');
    for i = 1:length(b)
        for j = 1:length(b{i}{1})
            bi = b{i}{1}{j};
            pi1 = sort(pp{i}{1}{bi(1,1)});
            pi2 = sort(pp{i}{1}{bi(end,1)}); 
            p.pattern{j}.occurrence{1}.start = ...
                fp{i}{1}(1,bi(1,1)) - mean(fp{i}{1}(1:2,pi1(bi(1,2))));
            p.pattern{j}.occurrence{2}.start = ...
                fp{i}{1}(1,bi(1,1));
            p.pattern{j}.occurrence{1}.end = ...
                fp{i}{1}(1,bi(end,1)) - mean(fp{i}{1}(1:2,pi2(bi(end,2))));
            p.pattern{j}.occurrence{2}.end = ...
                fp{i}{1}(1,bi(end,1));
        end
    end
    p = class(p,'mirpattern',mirdata(orig));
end
if postoption.period
    for i = 1:length(p.pattern)
        poi = p.pattern{i}.occurrence;
        if poi{1}.end > poi{2}.start
            poi{1}.end = poi{2}.start;
            cycle = poi{1}.end - poi{1}.start;
            ncycles = floor((poi{2}.end-poi{2}.start)/cycle)+2;
            poi{ncycles}.end = poi{2}.end;
            poi{2}.end = poi{2}.start + cycle;
            for j = 2:ncycles-1
                poi{j}.end = poi{j}.start + cycle;
                poi{j+1}.start = poi{j}.end;
            end
        end
        p.pattern{i}.occurrence = poi;
    end
end