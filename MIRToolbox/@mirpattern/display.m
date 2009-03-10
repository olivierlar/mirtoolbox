function display(p)

% MIRPATTERN/DISPLAY display of pattern

figure
hold on

k = 0;

for i = 1:length(p.pattern)
    display(['Pattern # ',num2str(i)])
    col = num2col(i);
    for j = 1:length(p.pattern{i}.occurrence)
        display(['Occurrence # ',num2str(j)])
        p.pattern{i}.occurrence{j}
        fill([p.pattern{i}.occurrence{j}.start ...
              p.pattern{i}.occurrence{j}.end ...
              p.pattern{i}.occurrence{j}.end ...
              p.pattern{i}.occurrence{j}.start], ...
              [k k k+1 k+1],col);
        k = k+1;
    end
    display('**************')
end

xlabel('time (in s.)')
ylabel('pattern and their occurrences')
set(gca,'YTick',[])