function varargout = mirpartial(orig,varargin)
             

        max.key = 'Max';
        max.type = 'Integer';
        max.default = Inf;
    option.max = max;
        
specif.option = option;

varargout = mirfunction(@mirpartial,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
type = 'mirpartial';


function p = main(orig,option,postoption)
a = get(orig,'Amplitude');
p = class(struct,'mirpartial',mirdata(orig));
p = purgedata(p);
%fp = get(orig,'FramePos');
pos = cell(1,length(a));
dat = cell(1,length(a));
tp = cell(1,length(a));
for i = 1:length(a)
    pos{i} = cell(1,length(a{i}));
    dat{i} = cell(1,length(a{i}));
    for j = 1:length(a{i})
        sizj = min(size(a{i}{j}{1},1),option.max);
        dat{i}{j} = a{i}{j}{1}(1:sizj,:);
        pos{i}{j} = repmat((1:sizj)',[1 size(dat{i}{j},2)]);
    end
end
p = set(p,'Title','Energy along partials',...
          'Abs','partials','Ord','energy','Data',dat,'Pos',pos,...
          'TrackPos',tp,'TrackVal',tp);