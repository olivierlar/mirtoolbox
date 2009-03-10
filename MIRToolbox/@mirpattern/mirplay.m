function varargout = mirplay(p,varargin)

        pat.key = 'Pattern';
        pat.type = 'Integer';
        pat.default = 0;
    option.pat = pat;
            
specif.option = option;

specif.eachchunk = 'Normal';

varargout = mirfunction(@mirplay,p,varargin,nargout,specif,@init,@main);
if nargout == 0
    varargout = {};
end


function [x type] = init(x,option)
type = '';


function noargout = main(p,option,postoption)
if not(option.pat)
    option.pat = 1:length(p.pattern);
end
n = get(p,'Name');
for h = 1:length(n)
    for i = option.pat
        display(['Pattern # ',num2str(i)])
        for j = 1:length(p.pattern{i}.occurrence)
            display(['Occurrence # ',num2str(j)])
            a = miraudio(n{h},'Extract',p.pattern{i}.occurrence{j}.start,...
                                     p.pattern{i}.occurrence{j}.end);
            mirplay(a)
        end
    end
end
noargout = {};