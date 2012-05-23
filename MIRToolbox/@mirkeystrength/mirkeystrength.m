function varargout = mirkeystrength(orig,varargin)
%   ks = mirkeystrength(x) computes the key strength, i.e., the probability
%   	associated with each possible key candidate.
%   Optional parameters:
%       mirkeystrength(...,'Frame',l,h) orders a frame decomposition of window
%           length l (in seconds) and hop factor h, expressed relatively to
%           the window length. For instance h = 1 indicates no overlap.
%           Default values: l = 1 seconds and h = .5
%       The mirchromagram options 'Weight' and 'Triangle' can be specified.
%   [ks,c] = mirkeystrength(...) also displays the chromagram used for the key 
%       strength estimation.
%
% Krumhansl, Cognitive foundations of musical pitch. Oxford UP, 1990.
% Gomez, Tonal description of polyphonic audio for music content processing,
%   INFORMS Journal on Computing, 18-3, pp. 294-304, 2006.

        wth.key = 'Weight';
        wth.type = 'Integer';
        wth.default = .5;
    option.wth = wth;
    
        tri.key = 'Triangle';
        tri.type = 'Boolean';
        tri.default = 0;
    option.tri = tri;
    
        transp.key = 'Transpose';
        transp.type = 'Integer';
        transp.default = 0;
        transp.when = 'After';
    option.transp = transp;

specif.option = option;
specif.defaultframelength = .1;
specif.defaultframehop = .125;

varargout = mirfunction(@mirkeystrength,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if not(isamir(x,'mirkeystrength'))
    if not(isamir(x,'mirchromagram'))
        x = mirchromagram(x,'Weight',option.wth,'Triangle',option.tri,'Normal');
    else
        x = mirchromagram(x,'Wrap','Normal');
    end
end
type = 'mirkeystrength';


function k = main(orig,option,postoption)
if iscell(orig)
    orig = orig{1};
end
if isa(orig,'mirkeystrength')
    c = [];
    k = orig;
else
    c = orig;
    load gomezprofs;
    m = get(c,'Magnitude');
    st = cell(1,length(m));
    kk = cell(1,length(m));
    %disp('Computing key strengths...')
    for i = 1:length(m)
        mi = m{i};
        if not(iscell(mi))
            mi = {mi};
        end
        si = cell(1,length(mi));
        ki = cell(1,length(mi));
        for j = 1:length(mi)
            mj = mi{j};
            sj = zeros(12,size(mj,2),size(mj,3),2);
            kj = cell(12,size(mj,2),size(mj,3));
            for k = 1:size(mj,2)
                for l = 1:size(mj,3)
                    if ~max(abs(mj(:,k,l)))
                        sj(:,k,l,:) = 0;
                    else
                        tmp = corrcoef([mj(:,k,l) gomezprofs']);
                        sj(:,k,l,1) = tmp(1,2:13);
                        sj(:,k,l,2) = tmp(1,14:25);
                    end
                    kj(:,k,l) = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
                end
            end
            si{j} = sj;
            ki{j} = kj;
        end
        st{i} = si;
        kk{i} = ki;
    end
    k = class(struct,'mirkeystrength',mirdata(c));
    k = purgedata(k);
    k = set(k,'Title','Key strength','Abs','tonal center','Ord','strength',...
              'Tonic',kk,'Strength',st,'MultiData',{'maj','min'},'Interpolable',0);
end
k = after(k,postoption);
k = {k c};


function k = after(k,postoption)
if postoption.transp
    transp = mod(postoption.transp,12);
    k = purgedata(k);
    d = mirgetdata(k);
    d = [d(13-transp:end,:,:,:);d(1:12-transp,:,:,:)];
    k = set(k,'Data',{{d}});
end