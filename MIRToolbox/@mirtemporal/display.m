function display(d,axis,songs)
% MIRDATA/DISPLAY display of a MIRtemporal

ST = dbstack;
if strcmp(ST(end).file,'arrayviewfunc.m')
    disp('To display its content in a figure, evaluate this variable directly in the Command Window.');
    return
end

if nargin<2
    axis = [];
end

if nargin<3
    songs = [];
end

if d.centered
    d = set(d,'Title',[get(d,'Title'),' (centered)']);
end
mirdisplay(mirdata(d),inputname(1),axis,songs);