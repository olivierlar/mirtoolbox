function display(d)
% MIRDATA/DISPLAY display of a MIRenvelope

if d.hwr
    d = set(d,'Title',[get(d,'Title'),' (half-wave rectified)']);
end
mirdisplay(mirtemporal(d),inputname(1));