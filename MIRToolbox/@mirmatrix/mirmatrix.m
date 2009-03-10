function m = mirmatrix(orig,varargin)
% Simple numerical matrix object used to display for instance correlation
% coefficients.

m = class(struct,'mirmatrix',mirdata(orig));
m = purgedata(m);
m = set(m,'FramePos',[],varargin{:});