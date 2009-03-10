function mirsave(d,varargin)
%   mirsave(d) saves data d in a file.
%   If d is a miraudio object,
%       misave(d) saves d as audio file(s). The file(s) name is based on 
%           the original file name(s), adding '.mir' before the standard
%           extension of the file.
%       mirsave(d,f) specifies the file names.
%           If d contains one single audio sequence, d is saved in a file
%               named f.
%           If d contains multiple audio sequences, each sequence is saved 
%               in a file whose name is the concatenation of the original
%               name and f.
%           If f ends with '.wav', the file is saved in WAV format (by
%               default).
%           If f ends with '.au', the file is saved in AU format.
%   If d is a mirscalar object,
%       mirsave(d) saves d in a text file, featuring three columns of data.
%           - the starting date of each frame (in seconds)
%           - the ending data of each frame (in seconds)
%           - the data itself.

mirsave(miraudio(d),varargin{:})