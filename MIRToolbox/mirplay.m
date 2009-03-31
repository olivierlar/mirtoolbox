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


        ch.key = 'Channel';
        ch.type = 'Integer';
        ch.default = 0;
    option.ch = ch;
        
        sg.key = 'Segment';
        sg.type = 'Integer';
        sg.default = 0;
    option.sg = sg;
        
        se.key = 'Sequence';
        se.type = 'Integer';
        se.default = 0;
    option.se = se;

        inc.key = 'Increasing';
        inc.type = 'MIRtb';
    option.inc = inc;

        dec.key = 'Decreasing';
        dec.type = 'MIRtb';
    option.dec = dec;

        every.key = 'Every';
        every.type = 'Integer';
    option.every = every;

specif.option = option;

specif.eachchunk = 'Normal';

varargout = mirfunction(@mirplay,a,varargin,nargout,specif,@init,@main);
if nargout == 0
    varargout = {};
end


function [x type] = init(x,option)
type = '';


function p = main(a,option,postoption)
mirplay(a);
p = {};