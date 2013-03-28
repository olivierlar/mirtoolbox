function varargout = maxcell(x,varargin)

specif.combinechunk = 'Average';

varargout = mirfunction(@maxcell,x,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
type = mirtype(a{1});


function y = main(x,option,postoption)
y = max(x{1},x{2});