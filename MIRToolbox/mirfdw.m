function res = mirfdw(filename,method,order)
%   res = mirfdtf(method,filename,order) computes spectral sub-band
% features using Filter Dependent Windowing (FDW).
%
% Input:
% - filename is the name of the audio file used as input
% - method is the method applied to each channel after the filterbank
% decomposition. By default, method = @mirflux
% - order is the filter order used as option in mirfilterbank(..., 'Manual'). 
% Default value: order = 4
%
% Initial algorithm developed by Fabi Prezja, University of Jyvaskyla, 2018
% Further adapted by Olivier Lartillot, 2019

    if nargin < 2
        method = @mirflux;
    end        
    if nargin < 3
        order = 4;
    end
    
    freq1 = [0,50,100,200,400,800,1600,3200,6400,12800,22050];
    freq2 = [-Inf,50,100,200,400,800,1600,3200,6400,12800,22050];
    for i = 1:length(freq1) - 1 
        sbi = mirfilterbank('Design', 'Manual', [freq2(i),freq2(i+1)], 'Order', order);
        fr = mirframe(sbi, frameconfig(freq1(i),freq1(i+1)), 50, '%');
        des = method(fr);
        resi = mireval(des, filename);
        res{i} = resi{1};
    end
end


function n = frameconfig(f1,f2)
%Determine the windowsize based on central frequency
%f1 = is the start frequency of the bandpass
%f2 = is the last frequency of  the bandapass
    if f2*f1==0, n=100/(f2*0.5);
    elseif f2 == +inf, print('Check f2')
    elseif f1 == -inf, print('Check f1')
    elseif f2/f1 >= 1.1, n=100/(sqrt(f2*f1));
    elseif f2/f1 < 1.1 , n=100/0.5*(f2+f1);
    else print('Something is wrong')
    end
end