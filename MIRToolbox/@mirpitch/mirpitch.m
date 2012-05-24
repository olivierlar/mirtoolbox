function varargout = mirpitch(orig,varargin)
%   p = mirpitch(x) evaluates the pitch frequencies (in Hz).
%   Specification of the method(s) for pitch estimation (these methods can
%       be combined):
%       mirpitch(...,'Autocor') computes an autocorrelation function
%           (Default method)
%           mirpitch(...'Enhanced',a) computes enhanced autocorrelation
%               (see help mirautocor)
%              toggled on by default
%           mirpitch(...,'Compress',k) performs magnitude compression
%               (see help mirautocor)
%           mirpitch(...,fb) specifies a type of filterbank.
%               Possible values:
%                   fb = 'NoFilterBank': no filterbank decomposition
%                   fb = '2Channels' (default value)
%                   fb = 'Gammatone' 
%       mirpitch(...,'AutocorSpectrum') computes the autocorrelation of
%           the FFT spectrum
%       mirpitch(...,'Cepstrum') computes the cepstrum
%       Alternatively, an autocorrelation or a cepstrum can be directly
%           given as first argument of the mirpitch function.
%   Peak picking options:
%       mirpitch(...,'Total',m) selects the m best pitches.
%           Default value: m = Inf, no limit is set concerning the number
%           of pitches to be detected.
%       mirpitch(...,'Mono') corresponds to morpitch(...,'Total',1)
%       mirpitch(...,'Min',mi) indicates the lowest frequency taken into
%           consideration.
%           Default value: 75 Hz. (Praat)
%       mirpitch(...,'Max',ma) indicates the highest frequency taken into
%           consideration. 
%           Default value: 2400 Hz. Because there seems to be some problems
%           with higher frequency, due probably to the absence of 
%           pre-whitening in our implementation of Tolonen and Karjalainen
%           approach (used by default, cf. below).
%       mirpitch(...,'Contrast',thr) specifies a threshold value.
%           (see help peaks)
%           Default value: thr = .1
%       mirpitch(...,'Order',o) specifies the ordering for the peak picking.
%           Default value: o = 'Amplitude'.
%       Alternatively, the result of a mirpeaks computation can be directly
%           given as first argument of the mirpitch function.
%   Post-processing options:
%       mirpitch(...,'Sum','no') does not sum back the channels at the end 
%           of the computation. The resulting pitch information remains
%           therefore decomposed into several channels.
%       mirpitch(...,'Median') performs a median filtering of the pitch
%           curve. When several pitches are extracted in each frame, the
%           pitch curve contains the best peak of each successive frame.
%       mirpitch(...,'Stable',th,n) remove pitch values when the difference 
%           (or more precisely absolute logarithmic quotient) with the
%           n precedent frames exceeds the threshold th. 
%           if th is not specified, the default value .1 is used
%           if n is not specified, the default value 3 is used
%       mirpitch(...'Reso',r) removes peaks whose distance to one or
%           several higher peaks is lower than a given threshold.
%           Possible value for the threshold r:
%               'SemiTone': ratio between the two peak positions equal to
%                   2^(1/12)
%       mirpitch(...,'Frame',l,h) orders a frame decomposition of window
%           length l (in seconds) and hop factor h, expressed relatively to
%           the window length. For instance h = 1 indicates no overlap.
%           Default values: l = 46.4 ms and h = 10 ms (Tolonen and
%           Karjalainen, 2000)
%   Preset model:
%       mirpitch(...,'Tolonen') implements (part of) the model proposed in
%           (Tolonen & Karjalainen, 2000). It is equivalent to
%           mirpitch(...,'Enhanced',2:10,'Generalized',.67,'2Channels')
%   [p,a] = mirpitch(...) also displays the result of the method chosen for
%       pitch estimation, and shows in particular the peaks corresponding
%       to the pitch values.
%   p = mirpitch(f,a,<r>) creates a mirpitch object based on the frequencies
%       specified in f and the related amplitudes specified in a, using a
%       frame sampling rate of r Hz (set by default to 100 Hz).
%
%   T. Tolonen, M. Karjalainen, "A Computationally Efficient Multipitch 
%       Analysis Model", IEEE TRANSACTIONS ON SPEECH AND AUDIO PROCESSING,
%       VOL. 8, NO. 6, NOVEMBER 2000

        ac.key = 'Autocor';
        ac.type = 'Boolean';
        ac.default = 0;
    option.ac = ac;
    
            enh.key = 'Enhanced';
            enh.type = 'Integer';
            enh.default = 2:10;
        option.enh = enh;

            filtertype.type = 'String';
            filtertype.choice = {'NoFilterBank','2Channels','Gammatone'};
            filtertype.default = '2Channels';
        option.filtertype = filtertype;

            sum.key = 'Sum';
            sum.type = 'Boolean';
            sum.default = 1;
        option.sum = sum;

            gener.key = {'Generalized','Compress'};
            gener.type = 'Integer';
            gener.default = .5;
        option.gener = gener;

        as.key = 'AutocorSpectrum';
        as.type = 'Boolean';
        as.default = 0;
    option.as = as;
    
        s.key = 'Spectrum';
        s.type = 'Boolean';
        s.default = 0;
    option.s = s;
        
        ce.key = 'Cepstrum';
        ce.type = 'Boolean';
        ce.default = 0;
    option.ce = ce;
        
%% peak picking options

        m.key = 'Total';
        m.type = 'Integer';
        m.default = Inf;
    option.m = m;
    
        multi.key = 'Multi';
        multi.type = 'Boolean';
        multi.default = 0;
    option.multi = multi;

        mono.key = 'Mono';
        mono.type = 'Boolean';
        mono.default = 0;
    option.mono = mono;

        mi.key = 'Min';
        mi.type = 'Integer';
        mi.default = 75;
    option.mi = mi;
        
        ma.key = 'Max';
        ma.type = 'Integer';
        ma.default = 2400;
    option.ma = ma;
        
        cthr.key = 'Contrast';
        cthr.type = 'Integer';
        cthr.default = .1;
    option.cthr = cthr;

        thr.key = 'Threshold';
        thr.type = 'Integer';
        thr.default = .4;
    option.thr = thr;

        order.key = 'Order';
        order.type = 'String';
        order.choice = {'Amplitude','Abscissa'};
        order.default = 'Amplitude';
    option.order = order;    

        reso.key = 'Reso';
        reso.type = 'String';
        reso.choice = {0,'SemiTone'};
        reso.default = 0;
    option.reso = reso;
        
        track.key = 'Track';        % Not used yet
        track.type = 'Boolean';
        track.default = 0;
    option.track = track;

%% post-processing options
        
        cent.key = 'Cent';
        cent.type = 'Boolean';
        cent.default = 0;
    option.cent = cent;
    
        segm.key = 'Segment';
        segm.type = 'Boolean';
        segm.default = 0;
    option.segm = segm;

        ref.key = 'Ref';
        ref.type = 'Integer';
        ref.default = 0;
    option.ref = ref;

        stable.key = 'Stable';
        stable.type = 'Integer';
        stable.number = 2;
        stable.default = [Inf 0];
        stable.keydefault = [.1 3];
    option.stable = stable;
    
        median.key = 'Median';
        median.type = 'Integer';
        median.default = 0;
        median.keydefault = .1;
    option.median = median;
    
        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [NaN NaN];
    option.frame = frame;
    
%% preset model

        tolo.key = 'Tolonen';
        tolo.type = 'Boolean';
        tolo.default = 0;
    option.tolo = tolo;
    
specif.option = option;
specif.chunkframebefore = 1;

if isnumeric(orig)
    if nargin<3
        f = 100;
    else
        f = varargin{2};
    end
    fp = (0:size(orig,1)-1)/f;
    fp = [fp;fp+1/f];
    p.amplitude = {{varargin{1}'}};
    s = mirscalar([],'Data',{{orig'}},'Title','Pitch','Unit','Hz',...
                     'FramePos',{{fp}},'Sampling',f,'Name',{inputname(1)});
    p = class(p,'mirpitch',s);
    varargout = {p};
else
    varargout = mirfunction(@mirpitch,orig,varargin,nargout,specif,@init,@main);
end



function [y type] = init(orig,option)
if option.tolo
    option.enh = 2:10;
    option.gener = .67;
    option.filtertype = '2Channels';
end
if not(option.ac) && not(option.as) && not(option.ce) && not(option.s)
    option.ac = 1;
end
if isnan(option.frame.length.val)
    option.frame.length.val = .0464;
end
if isnan(option.frame.hop.val)
    option.frame.hop.val = .01;
    option.frame.hop.unit = 's';
end
if isamir(orig,'mirscalar') || haspeaks(orig)
    y = orig;
else
    if isamir(orig,'mirautocor')
        y = mirautocor(orig,'Min',option.mi,'Hz','Max',option.ma,'Hz','Freq');
    elseif isamir(orig,'mircepstrum')
        y = orig;
    elseif isamir(orig,'mirspectrum')
        if not(option.as) && not(option.ce) && not(option.s)
            option.ce = 1;
        end
        if option.as
            y = mirautocor(orig,...
                            'Min',option.mi,'Hz','Max',option.ma,'Hz');
        end
        if option.ce
            ce = mircepstrum(orig,'freq',...
                            'Min',option.mi,'Hz','Max',option.ma,'Hz');
            if option.as
                y = y*ce;
            else
                y = ce;
            end
        end
        if option.s
            y = orig;
        end
    else
        if option.ac
            x = orig;
            if not(strcmpi(option.filtertype,'NoFilterBank'))
                x = mirfilterbank(x,option.filtertype);
            end
            x = mirframenow(x,option);
            y = mirautocor(x,'Generalized',option.gener);%,...
                               % 'Min',option.mi,'Hz','Max',option.ma,'Hz');
            if option.sum
                y = mirsummary(y);
            end
            y = mirautocor(y,'Enhanced',option.enh,'Freq');
            y = mirautocor(y,'Min',option.mi,'Hz','Max',option.ma,'Hz');
        end
        if option.as || option.ce || option.s
            x = mirframenow(orig,option);
            s = mirspectrum(x,'Min',option.mi,'Max',option.ma);
            if option.as
                as = mirautocor(s,...
                                'Min',option.mi,'Hz','Max',option.ma,'Hz');
                if option.ac
                    y = y*as;
                else
                    y = as;
                end
            end
            if option.ce
                ce = mircepstrum(s,'freq',...
                                'Min',option.mi,'Hz','Max',option.ma,'Hz');
                if option.ac || option.as
                    y = y*ce;
                else
                    y = ce;
                end
            end
            if option.s
                if option.ac || option.as
                    y = y*s;
                else
                    y = s;
                end
            end
        end
    end
end
type = {'mirpitch',mirtype(y)};
    

function o = main(x,option,postoption)
if option.multi && option.m == 1
    option.m = Inf;
end
if option.mono && option.m == Inf
    option.m = 1;
end
if iscell(x)
    if length(x)>1
        x2 = get(x{2},'Data');
        f2 = get(x{2},'Pos');
    end
    x = x{1};
else
    x2 = [];
end
if not(isa(x,'mirpitch'))
    x = mirpeaks(x,'Total',option.m,'Track',option.track,...
                   'Contrast',option.cthr,'Threshold',option.thr,...
                   'Reso',option.reso,'NoBegin','NoEnd',...
                   'Order',option.order);
end
if isa(x,'mirscalar')
    pf = get(x,'Data');
else
    pf = get(x,'PeakPrecisePos');
    pa = get(x,'PeakPreciseVal');
end
fp = get(x,'FramePos');

if option.cent
    for i = 1:length(pf)
        for j = 1:length(pf{i})
            for k = 1:size(pf{i}{j},3)
                for l = 1:size(pf{i}{j},2)    
                    pf{i}{j}{1,l,k} = 1200*log2(pf{i}{j}{1,l,k});
                end
            end
        end
    end
end

if option.segm
    for i = 1:length(pf)
        for j = 1:length(pf{i})
            for k = 1:size(pf{i}{j},3)
                startp = [];
                meanp = [];
                endp = [];
                deg = [];
                buffer = [];
                breaks = [];
                ref = option.ref;
                for l = 2:size(pf{i}{j},2)-1
                    if isempty(pf{i}{j}{1,l,k})
                        % Segment interrupted by no-pitch
                        if ~isempty(buffer)
                            meanp(end+1) = mean(buffer);
                            endp(end+1) = l-1;
                            [deg(end+1) ref] = cent2deg(meanp(end),ref);
                        end
                        buffer = [];
                        breaks(end+1) = l;
                    elseif isempty(buffer)
                        if 1 % abs(pf{i}{j}{1,l+1,k}-pf{i}{j}{1,l,k}) < 30
                            % New segment starting
                            startp(end+1) = l;
                            buffer = pf{i}{j}{1,l,k};
                        end
                    elseif abs(pf{i}{j}{1,l,k}-mean(buffer)) > 65
                        % Segment interrupted by pitch gap
                        meanp(end+1) = mean(buffer);
                        endp(end+1) = l-1;
                        [deg(end+1) ref] = cent2deg(meanp(end),ref);
                        if abs(pf{i}{j}{1,l+1,k}-pf{i}{j}{1,l,k}) < 30
                            % New segment starting
                            startp(end+1) = l;
                            buffer = pf{i}{j}{1,l,k};
                        else
                            buffer = [];
                        end
                    else
                        if length(pf{i}{j}{1,l,k})>1
                            mirerror('mirpitch','''Segment'' option only for monodies (use also ''Mono'')');
                        end
                        buffer(end+1) = pf{i}{j}{1,l,k};
                    end
                end
                
                if length(startp) > length(meanp)
                    startp(end) = [];
                end
                
                l = 1;
                while l <= length(endp)
                    if ~isempty(intersect(startp(l)-(1:5),breaks)) && ...
                            ~isempty(intersect(endp(l)+(1:5),breaks))
                        minlength = 8;
                    else
                        minlength = 2;
                    end
                    if endp(l)-startp(l) >= minlength
                    % Segment sufficiently long
                        found = 0;
                        pointer = startp(l);
                        for h = l-1:-1:max(l-10,1)
                            if deg(l) == deg(h)
                                if isempty(x2) % old version...
                                    test = ~isempty(pointer) && ...
                                            startp(l)-endp(h)<12 && ...
                                           ~any(intersect(startp(h):pointer-1,breaks));
                                    if test
                                        pointer = startp(h);
                                        test = abs(meanp(l)-meanp(h)) < 40;
                                    else
                                        pointer = [];
                                    end
                                else
                                    if 1 %abs(meanp(l)-meanp(h)) < 40;
                                        minp = min(meanp(l),meanp(h));
                                        maxp = max(meanp(l),meanp(h));
                                        minp = find(1200*log2(f2{i}{j}(:,1)) > minp,1);
                                        maxp = find(1200*log2(f2{i}{j}(:,1)) > maxp,1);
                                        zone = x2{i}{j}(minp:maxp,endp(h):startp(l));
                                        test = all(max(zone)>max(zone(:,1))*.05);
                                    end
                                    if test && size(zone,2)>2
                                        test = max(zone(:,end)) < max(max(zone(:,2:end-1)));
                                    end
                                end
                                if 0 %test
                                    % Segment close in frequency with recent one
                                    startp(l) = [];
                                    meanp(h) = mean(meanp([h l]));
                                    meanp(l) = [];
                                    deg(l) = [];
                                    endp(h) = endp(l);
                                    endp(l) = [];
                                    found = 1;
                                end
                                break
                            end
                        end
                        if ~found
                            l = l+1;
                        end
                    % Other cases: Segment too short
                    elseif l>1 && ...
                            startp(l) == endp(l-1)+1 && ...
                            abs(meanp(l)-meanp(l-1)) < 50
                        % Segment fused with previous one
                        startp(l) = [];
                        meanp(l-1) = mean(meanp(l-1:l));
                        meanp(l) = [];
                        deg(l) = [];
                        endp(l-1) = [];
                    elseif l < length(meanp) && ...
                            startp(l+1) == endp(l)+1 && ...
                            abs(meanp(l+1)-meanp(l)) < 50
                        % Segment fused with next one
                        startp(l+1) = [];
                        meanp(l) = mean(meanp(l:l+1));
                        meanp(l+1) = [];
                        deg(l+1) = [];
                        endp(l) = [];
                    else
                        % Segment removed
                        startp(l) = [];
                        meanp(l) = [];
                        deg(l) = [];
                        endp(l) = [];
                    end
                end
                
                ps{i}{j}{k} = startp;
                pe{i}{j}{k} = endp;
                pm{i}{j}{k} = meanp;
                dg{i}{j}{k} = deg;
            end
        end
    end
elseif isa(x,'mirpitch')
    ps = get(x,'Start');
    pe = get(x,'End');
    pm = get(x,'Mean');
    dg = get(x,'Degrees');
else
    ps = {};
    pe = {};
    pm = {};
    dg = {};
end

if option.stable(1) < Inf
    for i = 1:length(pf)
        for j = 1:length(pf{i})
            for k = 1:size(pf{i}{j},3)
                for l = size(pf{i}{j},2):-1:option.stable(2)+1
                    for m = length(pf{i}{j}{1,l,k}):-1:1
                        found = 0;
                        for h = 1:option.stable(2)
                            for n = 1:length(pf{i}{j}{1,l-h,k})
                                if abs(log10(pf{i}{j}{1,l,k}(m) ...
                                            /pf{i}{j}{1,l-h,k}(n))) ...
                                       < option.stable(1)
                                    found = 1;
                                end
                            end
                        end
                        if not(found)
                            pf{i}{j}{1,l,k}(m) = [];
                        end
                    end
                    pf{i}{j}{1,1,k} = zeros(1,0);
                end
            end
        end
    end
end
if option.median
    sr = get(x,'Sampling');
    for i = 1:length(pf)
        for j = 1:length(pf{i})
            if size(fp{i}{j},2) > 1
                npf = zeros(size(pf{i}{j}));
                for k = 1:size(pf{i}{j},3)
                    for l = 1:size(pf{i}{j},2)
                        if isempty(pf{i}{j}{1,l,k})
                            npf(1,l,k) = NaN;
                        else
                            npf(1,l,k) = pf{i}{j}{1,l,k}(1);
                        end
                    end
                end
                pf{i}{j} = medfilt1(npf,...
                     round(option.median/(fp{i}{j}(1,2)-fp{i}{j}(1,1))));
            end
        end
    end
end
if isa(x,'mirscalar')
    p.amplitude = 0;
else
    p.amplitude = pa;
end
p.start = ps;
p.end = pe;
p.mean = pm;
p.degrees = dg;
s = mirscalar(x,'Data',pf,'Title','Pitch','Unit','Hz');
p = class(p,'mirpitch',s);
o = {p,x};


function [deg ref] = cent2deg(cent,ref)
deg = round((cent-ref)/100);
%ref = cent - deg*100