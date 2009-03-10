function varargout = mirhcdf(orig,varargin)
%   df = mirhcdf(x) calculates the Harmonic Change Detection Function
%       related to x.
%   [df tc] = mirhcdf(x) also returns the tonal centroid vectors.
%   [df tc ch] = mirhcdf(x) also returns the intermediate chromagrams.
%
% C. A. Harte and M. B. Sandler, Detecting harmonic change in musical
%   audio, in Proceedings of Audio and Music Computing for Multimedia
%   Workshop, Santa Barbara, CA, 2006. 

specif.defaultframelength = .743;
specif.defaultframehop = .1;
varargout = mirfunction(@mirhcdf,orig,varargin,nargout,specif,@init,@main);


function [df type] = init(orig,option)
if isframed(orig)
    tc = mirtonalcentroid(orig);
else
    tc = mirtonalcentroid(orig,'Frame');
end
df = mirflux(tc);
type = 'mirscalar';


function df = main(df,option,postoption)