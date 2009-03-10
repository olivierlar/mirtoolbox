function display(c)

% CLASSIFY/DISPLAY display of classification

disp('Classification results:')
c.classes

if isnan(c.correct)
    disp('No label has been associated to the test set. Correct classification rate cannot be computed.');
else
    disp(['Correct classification rate: ',num2str(c.correct)]);
end

%disp(['Number of observations: ',num2str(c.nbobs)])
%disp(['Number of free parameters: ',num2str(c.nbparam)])

%disp('Posterior probability:')
%c.post{:}

if 0
    figure
    hold on
    vt = c.training;
    lt = c.labtraining;
    va = c.test;
    la = c.labtest;
    for i = 1:size(vt,2)
        scatter3(vt(1,i),vt(2,i),vt(3,i),'k+','SizeData',2);
        text(vt(1,i),vt(2,i),vt(3,i),lt{i},'Color','k');
    end
    for i = 1:size(va,2)
        scatter3(va(1,:),va(2,:),va(3,:),'r+','SizeData',2);
        text(va(1,i),va(2,i),va(3,i),la{i},'Color','r');
    end
    xlabel('1')
    ylabel('2')
    zlabel('3')
    rotate3d
end