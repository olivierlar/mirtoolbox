function varargout = mirplay(a,varargin)
%   mirplay(a) plays audio signal, envelope, or pitches.
%       If a is an envelope, what is actually played is a white noise of
%           same envelope.
%       If a is a mirpitch object, pitches are played using sinusoids.
%   Optional arguments:
%       mirplay(...,'Channel',i) plays the channel(s) of rank(s) indicated by 
%           the array i.
%       mirplay(...,'Segment',k) plays the segment(s) of rank(s) indicated by 
%           the array k.
%       mirplay(...,'Sequence',l) plays the sequence(s) of rank(s) indicated
%           by the array l.
%       mirplay(...,'Increasing',d) plays the sequences in increasing order
%           of d, which could be either an array or a mirscalar data.
%       mirplay(...,'Decreasing',d) plays the sequences in decreasing order
%           of d, which could be either an array or a mirscalar data.
%       mirplay(...,'Every',s) plays every s sequence, where s is a number
%           indicating the step between sequences.
%       Example: mirplay(mirenvelope('Folder'),...
%                        'increasing', mirrms('Folder'),...
%                        'every',5)

if ischar(a)
    varargout = mirplay(miraudio(a),varargin{:});
else
    mirerror('mirplay','You cannot play this type of object.')
end