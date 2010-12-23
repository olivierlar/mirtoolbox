function display(m)
% MIRMIDI/DISPLAY display of MIDI representation

figure
disp(' ');
d = get(m,'Data');
n = get(m,'Name');
for i = 1:length(d)
    pianoroll(d{i});
    fig = get(0,'CurrentFigure');
    va = inputname(1);
    if isempty(va)
        va = 'ans';
    end
    disp([va,' is the MIDI representation related to file ',n{i},...
        ' is displayed in Figure ',num2str(fig),'.']);
end
disp(' ');