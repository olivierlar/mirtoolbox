function display(d)
% MIRDESIGN/DISPLAY display of a MIR design

va = inputname(1);
if isempty(va)
    va = 'ans';
end
disp(' ');
for f = 1:length(d.fields)
    disp(['Field ',d.fields{f},' has value:']);
    d.data{f}
end
disp(['Temporary field has value:']);
d.tmp
