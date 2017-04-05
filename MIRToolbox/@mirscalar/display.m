function display(s)
% SCALAR/DISPLAY display the values of a scalar object

ST = dbstack;
if strcmp(ST(end).file,'arrayviewfunc.m')
    disp('To display its content in a figure, evaluate this variable directly in the Command Window.');
    return
end

display_figure(s);