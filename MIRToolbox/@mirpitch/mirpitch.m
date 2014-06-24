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
%       mirpitch(...,'Spectrum') computes the FFT spectrum
%       mirpitch(...,'AutocorSpectrum') computes the autocorrelation of
%           the FFT spectrum
%       mirpitch(...,'Cepstrum') computes the cepstrum
%       Alternatively, an autocorrelation or a cepstrum can be directly
%           given as first argument of the mirpitch function.
%   Peak picking options:
%       mirpitch(...,'Total',m) selects the m best pitches.
%           Default value: m = Inf, no limit is set concerning the number
%           of pitches to be detected.
%       mirpitch(...,'Mono') corresponds to mirpitch(...,'Total',1)
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
%       mirpitch(..., 'Cent') convert the pitch axis from Hz to cent scale.
%           One octave corresponds to 1200 cents, so that 100 cents
%           correspond to a semitone in equal temperament.
%       mirpitch(..., 'Segment') segments the obtained monodic pitch curve
%           in cents as a succession of notes with stable frequencies.
%           Additional parameters available: 'SegMinLength', 'SegPitchGap',
%               'SegTimeGap'.
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

            res.key = 'Res';
            res.type = 'Integer';
            res.default = NaN;
        option.res = res;
  
            db.key = 'dB';
            db.type = 'Integer';
            db.default = 0;
            db.keydefault = Inf;
        option.db = db;
        
        ce.key = 'Cepstrum';
        ce.type = 'Boolean';
        ce.default = 0;
    option.ce = ce;

        comb.key = 'Comb';
        comb.type = 'Boolean';
        comb.default = 0;
    option.comb = comb;
    
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
        segm.when = 'Both';
        segm.default = 0;
    option.segm = segm;

            segmin.key = 'SegMinLength';
            segmin.type = 'Integer';
            segmin.when = 'Both';
            segmin.default = 2;
        option.segmin = segmin;
        
            segpitch.key = 'SegPitchGap';
            segpitch.type = 'Integer';
            segpitch.when = 'Both';
            segpitch.default = 45;
        option.segpitch = segpitch;        

            segtime.key = 'SegTimeGap';
            segtime.type = 'Integer';
            segtime.when = 'Both';
            segtime.default = 20;
        option.segtime = segtime;      
        
            octgap.key = 'OctaveGap';
            octgap.type = 'Boolean';
            octgap.when = 'Both';
            octgap.default = 0;
        option.octgap = octgap;

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

        harmonic.key = 'Harmonic';
        harmonic.type = 'Boolean';
        harmonic.default = 0;
    option.harmonic = harmonic;
    
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
elseif option.harmonic
    option.s = 1;
    option.frame.hop.val = .1;
    option.res = 1;
    option.db = Inf;
end
if not(option.ac) && not(option.as) && not(option.ce) && not(option.s)
    option.ac = 1;
end
if option.segm && option.frame.length.val==0
    option.frame.length.val = NaN;
    option.frame.hop.val = NaN;
end
if isnan(option.frame.length.val)
    option.frame.length.val = .0464;
end
if isnan(option.frame.hop.val)
    option.frame.hop.val = .01;
    option.frame.hop.unit = 's';
end
if isamir(orig,'mirmidi') || isamir(orig,'mirscalar') || haspeaks(orig)
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
            if option.comb
                y = mirspectrum(x,'Min',option.mi,'Max',2000,'Res',1);%,'Sum');
            elseif option.s
                s = mirspectrum(x,'Min',option.mi,'Max',option.ma,...
                                  'Res',option.res,'dB',option.db);
                if option.ac
                    y = y*s;
                else
                    y = s;
                end
            end
            if option.as || option.ce
                s = mirspectrum(x);
                if option.as
                    as = mirautocor(s,'Min',option.mi,'Hz',...
                                      'Max',option.ma,'Hz');
                    if option.ac || option.s
                        y = y*as;
                    else
                        y = as;
                    end
                end
                if option.ce
                    ce = mircepstrum(s,'freq','Min',option.mi,'Hz',...
                                              'Max',option.ma,'Hz');
                    if option.ac || option.s || option.as
                        y = y*ce;
                    else
                        y = ce;
                    end
                end
            end
        end
    end
end
type = {'mirpitch',mirtype(y)};
    

function o = main(x,option,postoption)
if option.comb == 2
    option.m = Inf;
    option.order = 'Abscissa';
elseif option.multi && option.m == 1
    option.m = Inf;
elseif (option.mono && option.m == Inf) %|| option.segm
    option.m = 1;
elseif option.harmonic
    option.cthr = .01;
    option.thr = .5;
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

if option.comb == 1
    d = get(x,'Data');
    pos = get(x,'Pos');
    cb = cell(1,length(d));
    for i = 1:length(d)
        cb{i} = cell(1,length(d{i}));
        for j = 1:length(d{i})
            cb{i}{j} = zeros(size(d{i}{j},1),...
                             size(d{i}{j},2),...
                             size(d{i}{j},3));
            dij = d{i}{j}/max(max(max(d{i}{j})));
            for h = 1:size(d{i}{j},1)
                ph = pos{i}{j}(h,1,1);
                ip = h;
                for k = 2:size(d{i}{j},1)
                    [unused mp] = min(abs(pos{i}{j}(ip(end)+1:end,1,1) ...
                                          - ph * k));
                    if isempty(mp)
                        break
                    end
                    ip(end+1) = ip(end) + mp;
                end
                if length(ip) == 1
                    break
                end
                cbh = sum(dij(ip,:,:));
                for k = 1:length(ip)
                    cbh = cbh .* ...
                        (.5 * (2 - ...
                               exp(-(max(dij(ip(1:k),:,:),[],1).^2 * 5000))));
                end
                cb{i}{j}(h,:,:) = cbh;
            end
            cb{i}{j}(h+1:end,:,:) = [];
            pos{i}{j}(h+1:end,:,:) = [];
        end
    end
    x = set(x,'Data',cb,'Pos',pos,'Title','Spectral Comb');
end

if isa(x,'mirpitch')
    pf = get(x,'Data');
    pa = get(x,'Amplitude');
    if option.m < Inf
        for i = 1:length(pf)
            for j = 1:length(pf{i})
                for h = 1:length(pf{i}{j})
                    pf{i}{j}{h} = pf{i}{j}{h}(1:option.m,:);
                    pa{i}{j}{h} = pa{i}{j}{h}(1:option.m,:);
                end
            end
        end
    end
else
    if not(isa(x,'mirpitch') || isa(x,'mirmidi'))
        x = mirpeaks(x,'Total',option.m,'Track',option.track,...
                       'Contrast',option.cthr,'Threshold',option.thr,...
                       'Reso',option.reso,'NoBegin','NoEnd',...
                       'Order',option.order,'Harmonic',option.harmonic);
    end
    if isa(x,'mirscalar')
        pf = get(x,'Data');
    elseif option.harmonic
        pf = get(x,'TrackPos');
        pa = get(x,'TrackVal');
    else
        pf = get(x,'PeakPrecisePos');
        pa = get(x,'PeakPreciseVal');
    end
end

fp = get(x,'FramePos');

punit = 'Hz';

if option.comb == 2
    pp = get(x,'PeakPos');
    pv = get(x,'PeakVal');
    pm = get(x,'PeakMode');
    f = get(x,'Pos');
    for i = 1:length(pf)
        for j = 1:length(pf{i})
            maxf = f{i}{j}(end,1);
            for h = 1:length(pf{i}{j})
                sco = zeros(length(pf{i}{j}{h}),1);
                for k = 1:length(pf{i}{j}{h})
                    fk = pf{i}{j}{h}(k);
                    if fk > option.ma
                        break
                    end
                    ws = zeros(round(maxf / fk) ,1);
                    %err = mod(pf{i}{j}{h}/fk,1);
                    %err = min(err,1-err);
                    ws(1) = pa{i}{j}{h}(k);
                    for l = k+1:length(pf{i}{j}{h})
                        r = round(pf{i}{j}{h}(l) / fk);
                        if r == 1
                            continue
                        end
                        err = mod(pf{i}{j}{h}(l) / fk ,1);
                        err = min(err,1-err);
                        ws(r) = max(ws(r),pa{i}{j}{h}(l).*exp(-err^2*50));
                    end
                    
                    sco(k) = sum(ws);
                    if length(ws)>3 && ws(3)<.5
                        sco(k) = sco(k)/2;
                    end
                    %if length(ws)>5 && ws(5)<.5
                    %    sco(k) = sco(k)/2;
                    %end
                        %/(1+length(find(ws(2:end-1)<.01)));
                    %sco(k) = sum(pa{i}{j}{h}.*exp(-err));
                end
                %pa{i}{j}{h} = sco;
                [unused b] = max(sco);
                pf{i}{j}{h} = pf{i}{j}{h}(b);
                pa{i}{j}{h} = pa{i}{j}{h}(b);
                pp{i}{j}{h} = pp{i}{j}{h}(b);
                pv{i}{j}{h} = pv{i}{j}{h}(b);
                pm{i}{j}{h} = pm{i}{j}{h}(b);
            end
        end
    end
    x = set(x,'PeakPrecisePos',pf,'PeakPreciseVal',pa,...
              'PeakPos',pp,'PeakVal',pv,'PeakMode',pm);
end

if (option.cent || option.segm) && ...
        (~isa(x,'mirpitch') || strcmp(get(x,'Unit'),'Hz'))
    punit = 'cents';
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
    scale = [];
    for i = 1:length(pf)
        for j = 1:length(pf{i})
            if size(pf{i}{j},2) == 1 && size(pf{i}{j}{1},2) > 1
                pfj = cell(1,size(pf{i}{j}{1},2));
                paj = cell(1,size(pa{i}{j}{1},2));
                for l = 1:size(pf{i}{j}{1},2)
                    if isnan(pf{i}{j}{1}(l))
                        pfj{l} = [];
                        paj{l} = 0;
                    else
                        pfj{l} = pf{i}{j}{1}(l);
                        paj{l} = pa{i}{j}{1}(l);
                    end
                end
                pf{i}{j} = pfj;
                pa{i}{j} = paj;
            end
            
            for k = 1:size(pf{i}{j},3)
                startp = [];
                meanp = [];
                endp = [];
                deg = [];
                stabl = [];
                buffer = [];
                breaks = [];
                currentp = [];
                maxp = 0;
                reson = [];
                attack = [];
                
                if ~isempty(pf{i}{j}{1,1,k})
                    pf{i}{j}{1,1,k} = pf{i}{j}{1,1,k}(1);
                    pa{i}{j}{1,1,k} = pa{i}{j}{1,1,k}(1);
                end
                if ~isempty(pf{i}{j}{1,end,k})
                    pf{i}{j}{1,end,k} = pf{i}{j}{1,end,k}(1);
                    pa{i}{j}{1,end,k} = pa{i}{j}{1,end,k}(1);
                end
                
                for l = 2:size(pf{i}{j},2)-1
                    if ~isempty(pa{i}{j}{1,l,k}) && ...
                            pa{i}{j}{1,l,k}(1) > maxp
                        maxp = pa{i}{j}{1,l,k}(1);
                    end
                    if ~isempty(reson) && l-reson(1).end>50
                        reson(1) = [];
                    end
                    
                    if ~isempty(pf{i}{j}{1,l,k})
                        if 1 %isempty(pf{i}{j}{1,l-1,k})
                            pf{i}{j}{1,l,k} = pf{i}{j}{1,l,k}(1);
                            pa{i}{j}{1,l,k} = pa{i}{j}{1,l,k}(1);                            
                        else
                            [dpf idx] = min(abs(pf{i}{j}{1,l,k} - ...
                                                pf{i}{j}{1,l-1,k}));
                            if idx > 1 && ...
                                    pa{i}{j}{1,l,k}(1) - pa{i}{j}{1,l,k}(idx) > .02
                                pf{i}{j}{1,l,k} = pf{i}{j}{1,l,k}(1);
                                pa{i}{j}{1,l,k} = pa{i}{j}{1,l,k}(1);
                            else
                                pf{i}{j}{1,l,k} = pf{i}{j}{1,l,k}(idx);
                                pa{i}{j}{1,l,k} = pa{i}{j}{1,l,k}(idx);
                            end
                        end
                    end
                    
                    interrupt = 0;
                    if l == size(pf{i}{j},2)-1 || ...
                            isempty(pf{i}{j}{1,l,k}) || ...
                            (~isempty(buffer) && ...
                             abs(pf{i}{j}{1,l,k} - pf{i}{j}{1,l-1,k})...
                                > option.segpitch) || ...
                            (~isempty(currentp) && ...
                             abs(pf{i}{j}{1,l,k} - currentp) > ...
                                option.segpitch)
                        interrupt = 1;
                    elseif (~isempty(pa{i}{j}{1,l-1,k}) && ...
                             pa{i}{j}{1,l,k} - pa{i}{j}{1,l-1,k} > .01)
                        interrupt = 2;
                    end
                    
                    if ~interrupt
                        for h = 1:length(reson)
                            if abs(pf{i}{j}{1,l,k}-reson(h).pitch) < 50 && ...
                                    pa{i}{j}{1,l,k} < reson(h).amp/5
                                pa{i}{j}{1,l,k} = [];
                                pf{i}{j}{1,l,k} = [];
                                interrupt = 1;
                                break
                            end
                        end
                    end
                    
                    if interrupt
                        % Segment interrupted
                        if isempty(buffer) || ...
                                ...%length(buffer.pitch) < option.segmin || ...
                                0 %std(buffer.pitch) > 25
                             if length(startp) > length(endp)
                                startp(end) = [];
                            end
                        else
                            if isempty(currentp)
                                strong = find(buffer.amp > max(buffer.amp)*.75);
                                meanp(end+1) = mean(buffer.pitch(strong));
                            else
                                meanp(end+1) = currentp;
                            end
                            endp(end+1) = l-1;
                            hp = hist(buffer.pitch,5);
                            hp = hp/sum(hp);
                            entrp = -sum(hp.*log(hp+1e-12))./log(length(hp));
                            stabl(end+1) = entrp>.7;                        
                            deg(end+1) = cent2deg(meanp(end),scale);
                            reson(end+1).pitch = meanp(end);
                            reson(end).amp = mean(buffer.amp);
                            reson(end).end = l-1;
                            attack(end+1) = max(buffer.amp) > .05;
                        end
                              
                        if isempty(pf{i}{j}{1,l,k})
                            buffer = [];
                        else
                            buffer.pitch = pf{i}{j}{1,l,k};
                            buffer.amp = pa{i}{j}{1,l,k};
                            startp(end+1) = l;
                        end
                        currentp = [];
                        breaks(end+1) = l;
                        
                    elseif isempty(buffer)
                        % New segment starting
                        startp(end+1) = l;
                        buffer.pitch = pf{i}{j}{1,l,k};
                        buffer.amp = pa{i}{j}{1,l,k};
                        
                    else
                        if length(pf{i}{j}{1,l,k})>1
                            mirerror('mirpitch','''Segment'' option only for monodies (use also ''Mono'')');
                        end
                        buffer.pitch(end+1) = pf{i}{j}{1,l,k};
                        buffer.amp(end+1) = pa{i}{j}{1,l,k};
                        if length(buffer.pitch) > 4 && ...
                                std(buffer.pitch(1:end)) < 5 && ...
                                buffer.amp(end) > max(buffer.amp)*.5
                            currentp = mean(buffer.pitch(1:end));
                        %else
                        %    l
                        end
                    end                    
                end
                
                if length(startp) > length(meanp)
                    startp(end) = [];
                end
                
                l = 1;
                while l <= length(endp)
                    if 1 %~isempty(intersect(startp(l)-(1:5),breaks)) && ...
                         %  ~isempty(intersect(endp(l)+(1:5),breaks))
                        if 1 %attack(l)
                            minlength = option.segmin;
                        else
                            minlength = 6;
                        end
                    else
                        minlength = 2;
                    end
                    if endp(l)-startp(l) > minlength
                    % Segment sufficiently long
                        if l>1 && ~attack(l) && ...
                           startp(l) <= endp(l-1)+option.segtime && ...
                            abs(meanp(l)-meanp(l-1)) < 50
                                % Segment fused with previous one
                                startp(l) = [];
                                %meanp(l-1) = mean(meanp(l-1:l));
                                meanp(l) = [];
                                deg(l-1) = cent2deg(meanp(l-1),scale);
                                deg(l) = [];
                                attack(l-1) = max(attack(l),attack(l-1));
                                attack(l) = [];
                                endp(l-1) = [];
                                found = 1;
                        else
                                l = l+1;
                        end
                    % Other cases: Segment too short
                    elseif l>1 && ...
                            startp(l) <= endp(l-1)+option.segtime && ...
                            abs(meanp(l)-meanp(l-1)) < 50
                        % Segment fused with previous one
                        startp(l) = [];
                        %meanp(l-1) = mean(meanp(l-1:l));
                        meanp(l) = [];
                        deg(l) = [];
                        attack(l-1) = max(attack(l),attack(l-1));
                        attack(l) = [];
                        endp(l-1) = [];
                    elseif 0 && l < length(meanp) && ...
                            startp(l+1) <= endp(l)+option.segtime && ...
                            abs(meanp(l+1)-meanp(l)) < 50
                        % Segment fused with next one
                        startp(l+1) = [];
                        meanp(l) = meanp(l+1); %mean(meanp(l:l+1));
                        meanp(l+1) = [];
                        deg(l) = deg(l+1);
                        deg(l+1) = [];
                        attack(l) = max(attack(l),attack(l+1));
                        attack(l+1) = [];
                        endp(l) = [];
                    else
                        % Segment removed
                        startp(l) = [];
                        meanp(l) = [];
                        deg(l) = [];
                        attack(l) = [];
                        endp(l) = [];
                    end
                end               
                               
                l = 1;
                while l <= length(endp)
                    if (max([pa{i}{j}{1,startp(l):endp(l),k}]) < maxp/20 ...
                                && isempty(pa{i}{j}{1,startp(l)-1,k}) ...
                                && isempty(pa{i}{j}{1,endp(l)+1,k})) ...
                            || endp(l) - startp(l) < option.segmin
                        % Segment removed
                        fusetest = endp(l) - startp(l) < option.segmin;                        
                        startp(l) = [];
                        meanp(l) = [];
                        deg(l) = [];
                        endp(l) = [];
                        stabl(l) = [];
                        attack(l) = [];
                        
                        if fusetest && ...
                                l > 1 && l <= length(meanp) && ...
                                abs(meanp(l-1)-meanp(l)) < 50
                            % Preceding segment fused with next one
                            startp(l) = [];
                            meanp(l-1) = meanp(l); %mean(meanp(l:l+1));
                            meanp(l) = [];
                            deg(l-1) = deg(l);
                            deg(l) = [];
                            attack(l-1) = max(attack(l),attack(l-1));
                            attack(l) = [];
                            endp(l-1) = [];
                        end
                    else
                        l = l+1;
                    end
                end

                if option.octgap
                    l = 2;
                    while l <= length(endp)
                        if abs(meanp(l-1) - meanp(l) - 1200) < 50
                            % Segment removed
                            startp(l) = [];
                            meanp(l-1) = meanp(l);
                            meanp(l) = [];
                            deg(l-1) = deg(l);
                            deg(l) = [];
                            attack(l) = [];
                            endp(l-1) = [];
                            stabl(l) = [];
                        elseif abs(meanp(l) - meanp(l-1) - 1200) < 50
                            % Segment removed
                            startp(l) = [];
                            meanp(l) = meanp(l-1);
                            meanp(l) = [];
                            deg(l) = deg(l-1);
                            deg(l) = [];
                            attack(l) = [];
                            endp(l-1) = [];
                            stabl(l) = [];
                        else
                            l = l+1;
                        end
                    end
                end
                
                ps{i}{j}{k} = startp;
                pe{i}{j}{k} = endp;
                pm{i}{j}{k} = meanp;
                stb{i}{j}{k} = stabl;
                dg = {}; %{i}{j}{k} = deg;
            end
        end
    end
elseif isa(x,'mirpitch')
    ps = get(x,'Start');
    pe = get(x,'End');
    pm = get(x,'Mean');
    dg = get(x,'Degrees');
    stb = get(x,'Stable');
elseif isa(x,'mirmidi')
    nm = get(x,'Data');
    for i = 1:length(nm)
        startp = nm{i}(:,1);
        endp = startp + nm{i}(:,2);
        fp{i} = [startp endp]';
        ps{i} = {{1:length(startp)}};
        pe{i} = {{1:length(endp)}};
        pm{i} = {{nm{i}(:,4)'-68}};
        dg{i} = pm{i};
        stb{i} = [];
        pf{i} = {NaN(size(startp'))};
    end
    x = set(x,'FramePos',{fp});
else
    ps = {};
    pe = {};
    pm = {};
    dg = {};
    stb = {};
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
if 0 %isa(x,'mirscalar')
    p.amplitude = 0;
else
    p.amplitude = pa;
end
p.start = ps;
p.end = pe;
p.mean = pm;
p.degrees = dg;
p.stable = stb;
s = mirscalar(x,'Data',pf,'Title','Pitch','Unit',punit);
p = class(p,'mirpitch',s);
o = {p,x};


function [deg ref] = cent2deg(cent,ref)
deg = round((cent-ref)/100);
if isempty(deg)
    deg = 0;
end
%ref = cent - deg*100