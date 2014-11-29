function display(s)
% SCALAR/DISPLAY display the values of an histogram
disp(' ');
d = get(s,'Data');
b = get(s,'Bins');
t = ['Histogram of ',get(s,'Title')];
n = get(s,'Name');
l = get(s,'MultiData');
u = get(s,'Unit');
if not(isempty(u)) && not(strcmp(u,'.'))
    u = [' (in ',u,')'];
end
pp = get(s,'PeakPos');
for i = 1:length(d)
    figure
    set(gca,'NextPlot','replacechildren',...
            'LineStyleOrder',{'-',':','--',':','-.'})
    nl = size(b{i},1);
    for j = 1:nl
        bar((b{i}(j,:,1)+b{i}(j,:,2))/2,d{i}(j,:));
        hold all
    end
    if not(isempty(pp{i}{1}))
        for j = 1:nl;
            ppj = pp{i}{1}{1,j,1};
            bar((b{i}(j,ppj,1)+b{i}(j,ppj,2))/2,d{i}(j,ppj),'or')
        end
    end
    xlabel(['values',u]);
    ylabel('number of occurrences')
    title(t)
    axis tight
    axis 'auto y'
    nl = size(d{i},1);
    if nl>1
        legend(l,'Location','Best')
    end
    fig = get(0,'CurrentFigure');
    if isa(fig,'matlab.ui.Figure')
        fig = fig.Number;
    end
    va = inputname(1);
    if isempty(va)
        va = 'ans';
    end
    disp([va,' is the ',t,' related to file ',n{i},...
        ' is displayed in Figure ',num2str(fig),'.']);
end
disp(' ');