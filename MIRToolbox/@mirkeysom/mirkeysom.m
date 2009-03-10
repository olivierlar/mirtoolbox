function varargout = mirkeysom(orig,varargin)
%   ks = mirkeysom(x) projects a chromagram on a self-organizing map.
% Creates a pseudocolor map showing the projection of chromagram onto a 
% self-organizing map trained with the Krumhansl-Kessler profiles (modified
% for chromagrams). Colors correspond to Pearson correlation values.
% References:
%	Toiviainen, P. & Krumhansl, C. L. (2003). Measuring and modeling 
%	real-time responses to music: the dynamics of tonality induction. 
%	Perception, 32(6), 741-766.
%	Krumhansl, C. L., & Toiviainen, P. (2001) Tonal cognition. 
%	In R. J. Zatorre & I. Peretz (Eds.), The Biological Foundations of Music. 
%	Annals of the New York Academy of Sciences.

   %     filename.key = 'File';
   %     filename.type = 'String';
   %     filename.default = 0;
   % option.filename = filename;
    
specif.option = struct;
specif.defaultframelength = 1;
specif.defaultframehop = .5;

varargout = mirfunction(@mirkeysom,orig,varargin,nargout,specif,@init,@main);


function [c type] = init(orig,option)
c = mirchromagram(orig,'Normal');
type = 'mirkeysom';


function s = main(c,option,postoption)
if iscell(c)
    c = c{1};
end
load keysomaudiodata;
m = get(c,'Magnitude');
disp('Projecting the chromagram to a self-organizing map...')
z = cell(1,length(m));
p = cell(1,length(m));
for i = 1:length(m)
    mi = m{i};
    if not(iscell(mi))
        mi = {mi};
    end
    zi = cell(1,length(mi));
    pi = cell(1,length(mi));
    for j = 1:length(mi)
        mj = mi{j};
        zj = zeros(24,size(mj,2),size(mj,3),36);
        pi{j} = zeros(24,size(mj,2),size(mj,3),36);
        for k = 1:size(mj,2)
            for l = 1:size(mj,3)
                for kk=1:36
                    for ll=1:24
                        tmp = corrcoef([mj(:,k,l) somw(:,kk,ll)]);
                        zj(ll,k,l,kk) = tmp(1,2);
                    end
                end
            end
        end
        zi{j} = zj;
    end
    z{i} = zi;
    p{i} = pi;
end
s = class(struct,'mirkeysom',mirdata(c));
s = purgedata(s);
s = set(s,'Title','Key SOM','Abs','','Ord','','Weight',z,'Pos',p);
%if option.filename
%    mov = display(s,option.filename);
%else
%    mov = display(s);
%end
%s = {s,mov};