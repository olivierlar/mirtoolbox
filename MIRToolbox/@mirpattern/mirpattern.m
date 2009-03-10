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
    ap = get(orig,'AttackPos');
    rp = get(orig,'ReleasePos');
    fp = get(orig,'FramePos');
    p.pattern = {};
    for i = 1:length(ap)
        for j = 1:size(ap{i}{1},2)
            for k = 1:size(ap{i}{1},3)
                for l = 1:length(ap{i}{1}{1,j,k})
                    p.pattern{end+1}.occurrence{1}.start = ...
                        fp{i}{1}(1,ap{i}{1}{1,j,k}(l));
                    p.pattern{end}.occurrence{2}.start = ...
                        fp{i}{1}(1,ap{i}{1}{1,j,k}(l)) + mean(fp{i}{1}(1:2,j));
                    p.pattern{end}.occurrence{1}.end = ...
                        fp{i}{1}(2,rp{i}{1}{1,j,k}(l));
                    p.pattern{end}.occurrence{2}.end = ...
                        fp{i}{1}(2,rp{i}{1}{1,j,k}(l)) + mean(fp{i}{1}(1:2,j));
                end
            end
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