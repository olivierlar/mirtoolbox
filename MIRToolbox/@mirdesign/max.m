function varargout = max(a,b)

varargout = mirfunction(@maxcell,{a,b},{},1,struct,@init,@max);


function [x type] = init(x,option)
type = get(x{1},'Type');