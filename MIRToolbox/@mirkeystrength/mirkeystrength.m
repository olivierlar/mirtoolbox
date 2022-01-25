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
    
        meth.type = 'String';
        meth.choice = {'Gomez','Lartillot','Gomez+'};
        meth.default = 'Gomez';
    option.meth = meth;
    
        db.key = 'dB';
        db.type = 'Boolean';
        db.default = 0;
    option.db = db;
    
        origin.key = 'Tuning';
        origin.type = 'Integer';
        origin.default = 261.6256;
    option.origin = origin;
    
specif.option = option;
specif.defaultframelength = .1;
specif.defaultframehop = .125;

varargout = mirfunction(@mirkeystrength,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if not(isamir(x,'mirkeystrength'))
    if not(isamir(x,'mirchromagram'))
        x = mirchromagram(x,'Weight',option.wth,'Triangle',option.tri,'Normal',...
                            'dB',option.db,'Tuning',option.origin);
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
    if strcmpi(option.meth,'Gomez') || strcmpi(option.meth,'Gomez+')
        load gomezprofs;
        if strcmpi(option.meth,'Gomez+')
%             maj2 = [gomezprofs(1,1:8),gomezprofs(13,9:12)];
%             for i = 1:12
%                 gomezprofs(end+1,:) = maj2;
%                 maj2 = circshift(maj2,1);
%             end
            maj = [gomezprofs(13,[1,2,3,5,4,6,7,8,10,11,12,9])];
            for i = 1:12
                gomezprofs(end+1,:) = maj;
                maj = circshift(maj,1);
            end
            
        % figure,plot(0:11,gomezprofs(1,:))
        % set(gca,'Xtick',0:11)
        % set(gca,'XtickLabel',{'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'})
        % axis tight
        % figure,plot(0:11,gomezprofs(13,:))
        % set(gca,'Xtick',0:11)
        % set(gca,'XtickLabel',{'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'})
        % axis tight
        % figure,plot(0:11,gomezprofs(25,:))
        % set(gca,'Xtick',0:11)
        % set(gca,'XtickLabel',{'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'})
        % axis tight
            
        end
    else
        maj = [1 0 0 0 1 0 0 1 0 0 0 0]';
        min = [1 0 0 1 0 0 0 1 0 0 0 0]';
        majprofs = ones(12);
        minprofs = ones(12);
        for i = 1:12
            majprofs(:,i) = maj;
            minprofs(:,i) = min;
            maj = circshift(maj,1);
            min = circshift(min,1);
        end
    end
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
                    elseif strcmpi(option.meth,'Gomez') || strcmpi(option.meth,'Gomez+')
                        tmp = corrcoef([mj(:,k,l) gomezprofs']);
                        if strcmpi(option.meth,'Gomez+')
                            sj(:,k,l,1) = max(tmp(1,2:13),tmp(1,26:37));
                        else
                            sj(:,k,l,1) = tmp(1,2:13);
                        end
                        sj(:,k,l,2) = tmp(1,14:25);
                    else
                        sj(:,k,l,1) = mj(:,k,l)' * majprofs;
                        sj(:,k,l,2) = mj(:,k,l)' * minprofs;
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