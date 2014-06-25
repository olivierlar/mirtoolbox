function varargout = mirmfcc(orig,varargin)
%   c = mirmfcc(a) finds the Mel frequency cepstral coefficients (ceps),
%       a numerical description of the spectrum envelope.
%
%   Requires the Auditory Toolbox.
%
%   Optional arguments:
%       c = mirmfcc(...,'Rank',N) computes the coefficients of rank(s) N
%           (default: N = 1:13). Beware that the coefficient related to the
%           average energy is by convention here of rank 0. This zero value 
%           can be included to the array N as well.
%   If a is a frame decomposition, the temporal evolution of the MFCC,
%       along the successive frames, is returned. In this case, a second
%       option is available:
%       mirmfcc(...,'Delta',d) performs d temporal differentiations of
%           the coefficients, also called delta-MFCC (for d = 1) or
%           delta-delta-MFCC (for d = 2).
%       mirmfcc(...,'Delta') corresponds to mirmfcc(...,'Delta',1)
%   Optional arguments related to the delta computation:
%       mirmfcc(...,'Radius',r) specifies, for each frame, the number of
%           successive and previous neighbouring frames taken into
%           consideration for the least-square approximation.
%           Usually 1 or 2.
%           Default value: 2.

        nbbands.key = 'Bands';
        nbbands.type = 'Integer';
        nbbands.default = 40;
    option.nbbands = nbbands;

        rank.key = 'Rank';
        rank.type = 'Integer';
        rank.default = 1:13;
    option.rank = rank; 
        
        delta.key = 'Delta';
        delta.type = 'Integer';
        delta.default = 0;
        delta.keydefault = 1;
    option.delta = delta;
        
        radius.key = 'Radius';
        radius.type = 'Integer';
        radius.default = 2;
    option.radius = radius;
        
specif.option = option;
     
varargout = mirfunction(@mirmfcc,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if isamir(x,'miraudio') || isamir(x,'mirspectrum')
    x = mirspectrum(x,'Mel','log','Bands',option.nbbands);
end
type = 'mirmfcc';


function c = main(orig,option,postoption)
if iscell(orig)
    orig = orig{1};
end
if isa(orig,'mirmfcc')
    c = orig;
    if option.rank
        magn = get(c,'Data');
        rank = get(c,'Rank');
        for h = 1:length(magn)
            for k = 1:length(magn{h})
                m = magn{h}{k};
                r = rank{h}{k};
                r1 = r(:,1,1);
                range = find(ismember(r1,option.rank));
                magn{h}{k} = m(range,:,:);
                rank{h}{k} = r(range,:,:);
            end
        end
        c = set(c,'Data',magn,'Rank',rank);
    end
    c = modif(c,option);
else
    c.delta = 0;
    %disp('Computing Mel frequency cepstral coefficients...');
    e = get(orig,'Magnitude');

    % The following is largely based on the source code from Auditory Toolbox 
    % (A part that I could not call directly from MIRtoolbox)
    
    % (Malcolm Slaney, August 1993, (c) 1998 Interval Research Corporation)
    
    try
        MakeERBFilters(1,1,1); % Just to be sure that the Auditory Toolbox is installed
    catch
        error(['ERROR IN MIRFILTERBANK: Auditory Toolbox needs to be installed.']);
    end  
    
    dc = cell(1,length(e));
    rk = cell(1,length(e));
    for h = 1:length(e)
        dc{h} = cell(1,length(e{h}));
        rk{h} = cell(1,length(e{h}));
        for i = 1:length(e{h})
            ei = e{h}{i};
            totalFilters = size(ei,3); %Number of mel bands.

            % Figure out Discrete Cosine Transform.  We want a matrix
            % dct(i,j) which is totalFilters x cepstralCoefficients in size.
            % The i,j component is given by 
            %                cos( i * (j+0.5)/totalFilters pi )
            % where we have assumed that i and j start at 0.
            mfccDCTMatrix = 1/sqrt(totalFilters/2)*...
                            cos(option.rank' * ...
                                (2*(0:(totalFilters-1))+1) * ...
                                 pi/2/totalFilters);
            rank0 = find(option.rank == 0);
            mfccDCTMatrix(rank0,:) = mfccDCTMatrix(rank0,:) * sqrt(2)/2;
            ceps = zeros(size(mfccDCTMatrix,1),size(ei,2));
            for j = 1:size(ei,2)
                ceps(:,j) = mfccDCTMatrix * permute(ei(1,j,:),[3 1 2]);
            end
            dc{h}{i} = ceps;
            rk{h}{i} = repmat(option.rank(:),[1 size(ceps,2) size(ceps,3)]);
        end
    end
    c = class(c,'mirmfcc',mirdata(orig));
    c = purgedata(c);
    c = set(c,'Title','MFCC','Abs','coefficient ranks','Ord','magnitude',...
              'Data',dc,'Rank',rk);
    c = modif(c,option);
end
c = {c orig};


function c = modif(c,option)
d = get(c,'Data');
fp = get(c,'FramePos');
t = get(c,'Title');
if option.delta
    M = option.radius;
    for k = 1:option.delta
        for h = 1:length(d)
            for i = 1:length(d{h})
                nc = size(d{h}{i},2)-2*M;
                di = zeros(size(d{h}{i},1),nc);
                for j = 1:M
                    di = di + j * (d{h}{i}(:,M+j+(1:nc)) ...
                                 - d{h}{i}(:,M-j+(1:nc)));
                end
                di = di / 2 / sum((1:M).^2); % MULTIPLY BY 2 INSTEAD OF SQUARE FOR NORMALIZATION ?
                d{h}{i} = di;
                fp{h}{i} = fp{h}{i}(:,M+1:end-M);
            end
        end
        t = ['Delta-',t];
    end
end
c = set(c,'Data',d,'FramePos',fp,'Delta',get(c,'Delta')+option.delta,...
          'Title',t);