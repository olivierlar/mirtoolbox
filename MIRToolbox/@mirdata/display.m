function display(d,axis,songs,suffix)
% MIRDATA/DISPLAY display of a MIR data

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

if nargin<4
    suffix = [];
end

mirdisplay(d,inputname(1),axis,songs,suffix);