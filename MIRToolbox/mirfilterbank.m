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
        
%disp('Decomposing into a filterbank...');
d = get(x,'Data');
if size(d{1}{1},3) > 1
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
    if strcmpi(option.filtertype, 'Gammatone');
        [tmp x] = gettmp(x); %get(x,'Tmp');
        if not(Ch)
            Ch = 1:nCh;
        end
        erb = cell(1,length(d));
        for i = 1:length(d)
            erb{i} = cell(1,length(d{i}));
            nch{i} = Ch;
            for j = 1:length(d{i})
                try
                    coef = MakeERBFilters(f{i},nCh,option.lowF);
                catch
                    error(['ERROR IN MIRFILTERBANK: Auditory Toolbox needs to be installed.']);
                end                    
                [erb{i}{j} tmp] = ERBfilterbank(d{i}{j},coef,Ch,tmp);
            end
        end
        b = set(x,'Data',erb,'Channels',nch);%,'Tmp',tmp);
        clear erb
        b = settmp(b,tmp);
    elseif strcmpi(option.filtertype, '2Channels');
        if not(Ch)
            Ch = 1:2;
        end
        output = cell(1,length(d));
        try
            for i = 1:length(d)
                output{i} = cell(1,length(d{i}));
                [bl,al] = butter(4,[70 1000]/f{i}*2);
                k = 1;
                if ismember(1,Ch)
                    Hdl = dfilt.df2t(bl,al);
                    %fvtool(Hdl)
                    Hdl.PersistentMemory = true;
                    for j = 1:length(d{i})
                        low = filter(Hdl,d{i}{j});
                        output{i}{j}(:,:,k) = low;
                    end
                    k = k+1;
                end
                if ismember(2,Ch)
                    if f{i} < 20000
                        [bh,ah] = butter(2,1000/f{i}*2,'high');
                    else
                        [bh,ah] = butter(2,[1000 10000]/f{i}*2);
                    end
                    Hdlh = dfilt.df2t(bl,al);
                    %fvtool(Hdlh)
                    Hdh = dfilt.df2t(bh,ah);
                    %fvtool(Hdh)
                    Hdlh.PersistentMemory = true;
                    Hdh.PersistentMemory = true;
                    for j = 1:length(d{i})
                        high = filter(Hdh,d{i}{j});
                        high = max(high,0);
                        high = filter(Hdlh,high);
                        output{i}{j}(:,:,k) = high;
                    end
                end
                nch{i} = Ch;
            end
            b = set(x,'Data',output,'Channels',nch);
        catch
            warning(['WARNING IN MIRFILTERBANK: Signal Processing Toolbox (version 6.2.1, or higher) not installed: no filterbank decomposition.']);
            b = x;
        end
    elseif strcmpi(option.filtertype,'Manual');
        output = cell(1,length(d));
        for i = 1:length(d)
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
            output{i} = cell(1,length(d{i}));
            step = option.overlap;
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
                Hd(j) = dfilt.df2tsos(sos,g);
                Hd(j).PersistentMemory = true;
                for h = 1:length(d{i})
                    output{i}{h}(:,:,j) = filter(Hd(j),d{i}{h});
                end
            end
            %fvtool(Hd)
            nch{i} = 1:length(freqi)-step;
        end
        b = set(x,'Data',output,'Channels',nch);
    end
 end

    
function [output tmp] = ERBfilterbank(x, fcoefs, chans, tmp)
% The following is based on the source code from Auditory Toolbox 
% (A part that I could not call directly from MIRtoolbox, and that I
% generalize to chunk decomposition)

% (Malcolm Slaney, August 1993, (c) 1998 Interval Research Corporation)

% Process an input waveform with a gammatone filter bank.
% The fcoefs parameter, which completely specifies the Gammatone filterbank,
% should be designed with the MakeERBFilters function.

A0  = fcoefs(:,1);
A11 = fcoefs(:,2);
A12 = fcoefs(:,3);
A13 = fcoefs(:,4);
A14 = fcoefs(:,5);
A2  = fcoefs(:,6);
B0  = fcoefs(:,7);
B1  = fcoefs(:,8);
B2  = fcoefs(:,9);
gain= fcoefs(:,10);	
lc = length(chans);
output = zeros(size(x,1),size(x,2),lc);
if isempty(tmp)
    emptytmp = 1;
else
    emptytmp = 0;
end
for i = 1:lc
    chan = length(gain)-chans(i)+1; % Revert the channels order
    aa1 = [A0(chan)/gain(chan) A11(chan)/gain(chan)  A2(chan)/gain(chan)];
    bb1 = [B0(chan) B1(chan) B2(chan)];
    aa2 = [A0(chan) A12(chan) A2(chan)];
    bb2 = [B0(chan) B1(chan) B2(chan)];
    aa3 = [A0(chan) A12(chan) A2(chan)];
    bb3 = [B0(chan) B1(chan) B2(chan)];
    aa4 = [A0(chan) A13(chan) A2(chan)];
    bb4 = [B0(chan) B1(chan) B2(chan)];
    if emptytmp
        [y1 tmp(:,:,i,1)] = filter(aa1,bb1,x);
        [y2 tmp(:,:,i,2)] = filter(aa2,bb2,y1);
        [y3 tmp(:,:,i,3)] = filter(aa3,bb3,y2);
        [y4 tmp(:,:,i,4)] = filter(aa4,bb4,y3);
    else
        [y1 tmp(:,:,i,1)] = filter(aa1,bb1,x,tmp(:,:,i,1));
        [y2 tmp(:,:,i,2)] = filter(aa2,bb2,y1,tmp(:,:,i,2));
        [y3 tmp(:,:,i,3)] = filter(aa3,bb3,y2,tmp(:,:,i,3));
        [y4 tmp(:,:,i,4)] = filter(aa4,bb4,y3,tmp(:,:,i,4));
    end
    output(:,:,i) = y4;
end


function [y orig] = eachchunk(orig,option,missing)
y = mirfilterbank(orig,option);


function y = combinechunk(old,new)
do = get(old,'Data');
to = get(old,'Time');
dn = get(new,'Data');
tn = get(new,'Time');
y = set(old,'Data',{{[do{1}{1};dn{1}{1}]}},...
            'Time',{{[to{1}{1};tn{1}{1}]}});