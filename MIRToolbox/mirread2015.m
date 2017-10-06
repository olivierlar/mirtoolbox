function [d,tp,fp,f,l,b,n,ch] = mirread2015(extract,orig,load,folder,verbose)
% Read the audio file ORIG, at temporal position indicated by EXTRACT. If
% EXTRACT is empty, all the audio file is loaded.
%   If LOAD is set to 0, just the meta-data is collected, and the actual
%       audio data is not taken into consideration. If it is set to 1, the
%       data are loaded from the current directory. If LOAD is a string, it
%       is considered as the path of the directory.   
%   If FOLDER is set to 1, no error is returned if an audio file cannot be
%       loaded.
% Output:
%   D is the audio signal,
%   TP are the temporal positions,
%   FP are the two extreme temporal positions (used for frame positions),
%   F is the sampling rate,
%   L is the duration in seconds,
%   B is the resolution in number of bits,
%   N is the file name.
%   CH are the channel index.

if nargin < 5
    verbose = 0;
end
d = {};
f = {};
l = {};
b = {};
tp = {};
fp = {};
n = {};
ch = {};

if strcmp(orig,'.') || strcmp(orig,'..')
    return
end

try
    [d,f,l,b,tp,fp,n,ch] = audioreader(extract,orig,load,verbose,folder);
catch ME
    if folder
        return
    end
    rethrow(ME)
end

        
function [d,f,l,b,tp,fp,n,ch] = audioreader(extract,file,load,verbose,folder)
n = file;
if folder
    file = ['./',file];
end
if load
    if isempty(extract)
        [s,f] = audioread(file);
        i = audioinfo(file);
        if isfield(i,'BitsPerSample')
            b = i.BitsPerSample;
        elseif isfield(i,'BitRate')
            b = i.BitRate;
        else
            b = NaN;
        end
    else
        i = audioinfo(file);
        f = i.SampleRate;
        if isfield(i,'BitsPerSample')
            b = i.BitsPerSample;
        elseif isfield(i,'BitRate')
            b = i.BitRate;
        else
            b = NaN;
        end
        s = audioread(file,extract(1:2));
        if length(extract) > 2
            s = s(:,extract(3));
        end
    end
    if verbose
        disp([file,' loaded.']);
    end
    d{1} = reshape(s,size(s,1),1,size(s,2)); %channels along dim 3
    ch = 1:size(s,2);
    if isempty(extract)
        tp{1} = (0:size(s,1)-1)'/f;
    else
        tp{1} = (extract(1)-1+(0:size(s,1)-1))'/f;
    end
    l{1} = (size(s,1)-1)/f;
    if isempty(s)
        fp{1} = 0;
    else
        fp{1} = tp{1}([1 end]);
    end
else
    i = audioinfo(file);
    d = i.TotalSamples;
    f = i.SampleRate;
    if isfield(i,'BitsPerSample')
        b = i.BitsPerSample;
    elseif isfield(i,'BitRate')
        b = i.BitRate;
    else
        b = NaN;
    end
    ch = i.NumChannels;
    l = d/f;
    tp = {};
    fp = {};
end