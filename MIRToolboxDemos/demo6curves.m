function demo6curves(file)
% Example of use of curve analysis tools, such as moment estimations.

if nargin == 0
    file = 'trumpet.wav';
end

s = mirspectrum(file)
mircentroid(s)
mirspread(s)
mirskewness(s)
mirkurtosis(s)
mirflatness(s)
mirregularity(s)
mirentropy(s)
clear s

e = mirenvelope(file)
mircentroid(e)
mirspread(e)
mirskewness(e)
mirkurtosis(e)
mirflatness(e)
mirentropy(e)

display('Strike any key to continue...');
pause
close all

s = mirspectrum(file,'Frame',.1,1)
mircentroid(s)
mirspread(s)
mirskewness(s)
mirkurtosis(s)
mirflatness(s)
mirregularity(s)
mirentropy(s)
clear s

display('Strike any key to continue...');
pause
close all

fe = mirframe(e,.1,1)
clear e
mircentroid(fe)
mirspread(fe)
mirskewness(fe)
mirkurtosis(fe)
mirflatness(fe)
mirentropy(fe)
clear fe

display('Strike any key to continue...');
pause
close all

if nargin == 0
    file = 'czardas.wav';
end
t = mirtempo(file,'Periodicity','Frame')
h = mirhisto(t)
clear t
mircentroid(h)
mirspread(h)
mirskewness(h)
mirkurtosis(h)
mirflatness(h)