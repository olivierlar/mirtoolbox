function varargout = mirfilterbank(orig,varargin)
%   b = mirfilterbank(x) performs a filterbank decomposition of an audio
%       waveform.
%   Optional arguments:
%       mirfilterbank(...,t) selects a type of filterbank.
%           Possible values:
%               t = 'Gammatone' (default for audio files).
%                   Requires the Auditory Toolbox.
%                   mirfilterbank(...,'Lowest',f): lowest frequency in Hz
%                       (default: 50)
%               t = '2Channels' proposed in (Tolonen & Karjalainen, 2000)
%           mirfilterbank(...,'NbChannels',N), or simply filterbank(x,N):
%               specifies the number of channels in the bank.
%               (default: N = 10)
%           mirfilterbank(...,'Channel',c) only output the channels whose
%               ranks are indicated in the array c.
%               (default: c = (1:N))
%       mirfilterbank(...,'Manual',f) specifies a set of low-pass, band-
%                   pass and high-pass eliptic filters (Scheirer, 1998). 
%                   The series of cut-off frequencies as to be specified
%                   as next parameter f.
%                       If this series of frequencies f begins with -Inf,
%                           the first filter is low-pass.
%                       If this series of frequencies f ends with Inf,
%                           the last filter is high-pass.
%           mirfilterbank(...,'Order',o) specifies the order of the filter.
%               (Default: o = 4) (Scheirer, 1998)
%           mirfilterbank(...,'Hop',h) specifies the degree of spectral
%               overlapping between successive channels.
%               If h = 1 (default value), the filters are non-overlapping.
%               If h = 2, the filters are half-overlapping.
%               If h = 3, the spectral hop factor between successive
%                   filters is a third of the whole frequency region, etc.
%       mirfilterbank(...,p) specifies predefined filterbanks, all
%           implemented using elliptic filters, by default of order 4.
%           Possible values:
%               p = 'Mel' (mel-scale)
%               p = 'Bark' (bark-scale)
%               p = 'Scheirer' proposed in (Scheirer, 1998) corresponds to
%                   'Manual',[-Inf 200 400 800 1600 3200 Inf]
%               p = 'Klapuri' proposed in (Klapuri, 1999) corresponds to
%                   'Manual',44*[2.^ ([ 0:2, ( 9+(0:17) )/3 ]) ]

        nCh.key = 'NbChannels';
        nCh.type = 'Integer';
        nCh.default = 10;
    option.nCh = nCh;
    
        Ch.key = {'Channel','Channels'};
        Ch.type = 'Integer';
        Ch.default = 0;
    option.Ch = Ch;
    
        lowF.key = 'Lowest';
        lowF.type = 'Integer';
        lowF.default = 50;
    option.lowF = lowF;

        freq.key = 'Manual';
        freq.type = 'Integer';
        freq.default = NaN;
    option.freq = freq;
    
        overlap.key = 'Hop';
        overlap.type = 'Boolean';
        overlap.default = 1;
    option.overlap = overlap;

        filtertype.type = 'String';
        filtertype.choice = {'Gammatone','2Channels',0};
        filtertype.default = 'Gammatone';
    option.filtertype = filtertype;

        filterorder.key = 'Order';
        filterorder.type = 'Integer';
        filterorder.default = 4;
    option.filterorder = filterorder;

        presel.type = 'String';
        presel.choice = {'Scheirer','Klapuri','Mel','Bark'};
        presel.default = '';
    option.presel = presel;

specif.option = option;

specif.eachchunk = @eachchunk;
specif.combinechunk = @combinechunk;

varargout = mirfunction(@mirfilterbank,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
type = 'miraudio';


function b = main(x,option,postoption)
if iscell(x)
    x = x{1};
end

if ~isamir(x,'miraudio')
    mirerror('MIRFILTERBANK','The input should be an audio waveform.');
end

f = get(x,'Sampling');
    
if strcmpi(option.presel,'Scheirer')
    option.freq = [-Inf 200 400 800 1600 3200 Inf];
elseif strcmpi(option.presel,'Klapuri')
    option.freq = 44*[2.^([0:2,(9+(0:17))/3])];
elseif strcmpi(option.presel,'Mel')
    lowestFrequency = 133.3333;
    linearFilters = 13;
    linearSpacing = 66.66666666;
    logFilters = 27;
    logSpacing = 1.0711703;
    totalFilters = linearFilters + logFilters;
    cepstralCoefficients = 13;
    option.freq = lowestFrequency + (0:linearFilters-1)*linearSpacing;
    option.freq(linearFilters+1:totalFilters+2) = ...
        option.freq(linearFilters) * logSpacing.^(1:logFilters+2);

    option.overlap = 2;
elseif strcmpi(option.presel,'Bark')
    option.freq = [10 20 30 40 51 63 77 92 108 127 148 172 200 232 ...
                    270 315 370 440 530 640 770 950 1200 1550]*10; %% Hz
end
if not(isnan(option.freq))
    option.filtertype = 'Manual';
end
        
d = get(x,'Data');
if size(d{1}{1},3) > 1
    warning('WARNING IN MIRFILTERBANK: The input data is already decomposed into channels. No more channel decomposition.');
    if option.Ch
        cx = get(x,'Channels');
        db = cell(1,length(d));
        for k = 1:length(d)
            db{k} = cell(1,length(d{k}));
            for l = 1:length(d{k})
                for i = 1:length(option.Ch)
                    db{k}{l}(:,:,i) = d{k}{l}(:,:,find(cx{k} == option.Ch(i)));
                end
            end
        end
        b = set(x,'Data',db,'Channels',option.Ch);
    else
        b = x;
    end
else
    i = 1;
    while i <= length(d)
        if size(d{i}{1},2) > 1
            warning('WARNING IN MIRFILTERBANK: The frame decomposition should be performed after the filterbank decomposition.');
            disp(['Suggestion: Use the ' char(34) 'Frame' char(34) ' option available to some of the MIRtoolbox functions.'])
            break
        end
        i = i+1;
    end
    nCh = option.nCh;
    Ch = option.Ch;
    
    [tmp x] = gettmp(x);
    output = cell(1,length(d));
    nch = cell(1,length(d));
    for i = 1:length(d)
        if isempty(tmp)
            %% Determination of the filterbank specifications
            if strcmpi(option.filtertype,'Gammatone')
                if not(Ch)
                    Ch = 1:nCh;
                end
                Hd = ERBFilters(f{i},nCh,Ch,option.lowF);
                nch{i} = Ch;
            elseif strcmpi(option.filtertype,'2Channels')
                if not(Ch)
                    Ch = 1:2;
                end
                [bl,al] = butter(4,[70 1000]/f{i}*2);
                if ismember(1,Ch)
                    Hd{1} = dfilt.df2t(bl,al);
                    k = 2;
                else
                    k = 1;
                end
                if ismember(2,Ch)
                    if f{i} < 20000
                        [bh,ah] = butter(2,1000/f{i}*2,'high');
                    else
                        [bh,ah] = butter(2,[1000 10000]/f{i}*2);
                    end
                    Hd{k} = {dfilt.df2t(bl,al),...
                             @(x) max(x,0),...
                             dfilt.df2t(bh,ah)};
                end
                nch{i} = Ch;
            elseif strcmpi(option.filtertype,'Manual')
                freqi = option.freq;
                j = 1;
                while j <= length(freqi)
                    if not(isinf(freqi(j))) && freqi(j)>f{i}/2
                        if j == length(freqi)
                            freqi(j) = Inf;
                        else
                            freqi(j) = [];
                            j = j-1;
                        end
                    end
                    j = j+1;
                end
                step = option.overlap;
                if length(freqi) <= step
                    Hd = {};
                    nch = [];
                else
                    for j = 1:length(freqi)-step
                        if isinf(freqi(j))
                            [z{j},p{j},k{j}] = ellip(option.filterorder,3,40,...
                                                freqi(j+step)/f{i}*2);
                        elseif isinf(freqi(j+step))
                            [z{j},p{j},k{j}] = ellip(option.filterorder,3,40,...
                                                freqi(j)/f{i}*2,'high');
                        else
                            [z{j},p{j},k{j}] = ellip(option.filterorder,3,40,...
                                                freqi([j j+step])/f{i}*2);
                        end
                    end
                    for j = 1:length(z)
                        [sos,g] = zp2sos(z{j},p{j},k{j});
                        Hd{j} = dfilt.df2tsos(sos,g);
                    end
                    nch{i} = 1:length(freqi)-step;
                end
            end
            
            if length(d) == 1
                for k = 1:length(Hd)
                    Hdk = Hd{k};
                    if ~iscell(Hdk)
                        Hdk = {Hdk};
                    end
                    for h = 1:length(Hdk)
                        if ~isa(Hdk{h},'function_handle')
                            Hdk{h}.PersistentMemory = true;
                        end
                    end
                end
            end
        elseif i == 1
            Hd = tmp;
        end 
        
        output{i} = cell(1,length(d{i}));
        for j = 1:length(d{i})
            for k = 1:length(Hd)
                dk = d{i}{j};
                Hdk = Hd{k};
                if ~iscell(Hdk)
                    Hdk = {Hdk};
                end
                for h = 1:length(Hdk)
                    if isa(Hdk{h},'function_handle')
                        dk = Hdk{h}(dk);
                    else
                        dk = Hdk{h}.filter(dk);
                    end
                end
                output{i}{j}(:,:,k) = dk;
            end
        end
    end
    b = set(x,'Data',output,'Channels',nch);   
    b = settmp(b,Hd);
end

%%

function Hd=ERBFilters(fs,numChannels,chans,lowFreq)
% This function computes the filter coefficients for a bank of 
% Gammatone filters.  These filters were defined by Patterson and 
% Holdworth for simulating the cochlea.  
% The transfer function  of these four second order filters share the same
% denominator (poles) but have different numerators (zeros).
% The filter bank contains "numChannels" channels that extend from
% half the sampling rate (fs) to "lowFreq".

% Note this implementation fixes a problem in the original code by
% computing four separate second order filters.  This avoids a big
% problem with round off errors in cases of very small cfs (100Hz) and
% large sample rates (44kHz).  The problem is caused by roundoff error
% when a number of poles are combined, all very close to the unit
% circle.  Small errors in the eigth order coefficient, are multiplied
% when the eigth root is taken to give the pole location.  These small
% errors lead to poles outside the unit circle and instability.  Thanks
% to Julius Smith for leading me to the proper explanation.

% Code taken from Auditory Toolbox and optimized.
% (Malcolm Slaney, August 1993, (c) 1998 Interval Research Corporation)

T = 1/fs;
EarQ = 9.26449;				%  Glasberg and Moore Parameters
minBW = 24.7;

%%
% Computes an array of numChannels frequencies uniformly spaced between
% fs/2 and lowFreq on an ERB scale.
%
% For a definition of ERB, see Moore, B. C. J., and Glasberg, B. R. (1983).
% "Suggested formulae for calculating auditory-filter bandwidths and
% excitation patterns," J. Acoust. Soc. Am. 74, 750-753.
% 
% Derived from Apple TR #35, "An
% Efficient Implementation of the Patterson-Holdsworth Cochlear
% Filter Bank."  See pages 33-34.
cf = -(EarQ*minBW) + exp((1:numChannels)'*(-log(fs/2 + EarQ*minBW) + ...
		log(lowFreq + EarQ*minBW))/numChannels) * (fs/2 + EarQ*minBW);

%%
ERB = ((cf/EarQ) + minBW);
B=1.019*2*pi*ERB;

A0 = T;
A2 = 0;
B0 = 1;
B1 = -2*cos(2*cf*pi*T)./exp(B*T);
B2 = exp(-2*B*T);

A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;
A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;
A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;
A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;

gain = abs((-2*exp(4i*cf*pi*T)*T + ...
                 2*exp(-(B*T) + 2i*cf*pi*T).*T.* ...
                         (cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
                          sin(2*cf*pi*T))) .* ...
           (-2*exp(4i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * ...
               sin(2*cf*pi*T))).* ...
           (-2*exp(4i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) - ...
               sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
           (-2*exp(4i*cf*pi*T)*T + 2*exp(-(B*T) + 2i*cf*pi*T).*T.* ...
           (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
          (-2 ./ exp(2*B*T) - 2*exp(4i*cf*pi*T) +  ...
           2*(1 + exp(4i*cf*pi*T))./exp(B*T)).^4);
	
allfilts = ones(length(cf),1);
A0 = A0*allfilts;
A2 = A2*allfilts;
B0 = B0*allfilts;

for i = 1:length(chans)
    chan = length(gain)-chans(i)+1; % Revert the channels order
    aa1 = [A0(chan)/gain(chan) A11(chan)/gain(chan)  A2(chan)/gain(chan)];
    bb1 = [B0(chan) B1(chan) B2(chan)];
    aa2 = [A0(chan) A12(chan) A2(chan)];
    bb2 = [B0(chan) B1(chan) B2(chan)];
    aa3 = [A0(chan) A13(chan) A2(chan)];
    bb3 = [B0(chan) B1(chan) B2(chan)];
    aa4 = [A0(chan) A14(chan) A2(chan)];
    bb4 = [B0(chan) B1(chan) B2(chan)];
    Hd{i} = {dfilt.df2t(aa1,bb1);dfilt.df2t(aa2,bb2);...
             dfilt.df2t(aa3,bb3);dfilt.df2t(aa4,bb4)};
end


%%
function [y orig] = eachchunk(orig,option,missing)
y = mirfilterbank(orig,option);


function y = combinechunk(old,new)
do = get(old,'Data');
to = get(old,'Time');
dn = get(new,'Data');
tn = get(new,'Time');
y = set(old,'Data',{{[do{1}{1};dn{1}{1}]}},...
            'Time',{{[to{1}{1};tn{1}{1}]}});