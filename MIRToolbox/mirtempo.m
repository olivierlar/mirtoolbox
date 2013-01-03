function varargout = mirtempo(x,varargin)
%   t = mirtempo(x) evaluates the tempo in beats per minute (BPM).
%   Optional arguments:
%       mirtempo(...,'Total',m) selects not only the best tempo, but the m
%           best tempos.
%       mirtempo(...,'Frame',l,h) orders a frame decomposition of window
%           length l (in seconds) and hop factor h, expressed relatively to
%           the window length. For instance h = 1 indicates no overlap.
%           Default values: l = 3 seconds and h = .1
%       mirtempo(...,'Min',mi) indicates the lowest tempo taken into
%           consideration, expressed in bpm.
%           Default value: 40 bpm.
%       mirtempo(...,'Max',ma) indicates the highest tempo taken into
%           consideration, expressed in bpm.
%           Default value: 200 bpm.
%       mirtempo(...,s) selects the tempo estimation strategy:
%           s = 'Autocor': Approach based on the computation of the
%               autocorrelation. (Default strategy)
%               Option associated to the mirautocor function can be
%               passed here as well (see help mirautocor):
%                   'Enhanced' (toggled on by default here)
%           s = 'Spectrum': Approach based on the computation of the
%               spectrum .
%               Option associated to the mirspectrum function can be
%               passed here as well (see help mirspectrum):
%                   'ZeroPad' (set by default to 10000 samples)
%                   'Prod' (toggled off by default)
%           These two strategies can be combined: the autocorrelation
%               function is translated into the frequency domain in order
%               to be compared to the spectrum curve.
%               tempo(...,'Autocor','Spectrum') multiplies the two curves.
%           Alternatively, an autocorrelation function ac or a spectrum sp
%               can be directly passed to the function tempo:
%                   mirtempo(ac) or mirtempo(sp)
%       The options related to the onset detection phase can be specified 
%               here as well (see help mironsets):
%               onset detection strategies: 'Envelope', 'DiffEnvelope'
%               (corresponding to 'Envelope', 'Diff'), 'SpectralFlux,
%               'Pitch', 'Log', 'Mu', 'Filterbank'
%               mironsets(...,'Sum',w) specifies when to sum the channels.
%                   Possible values:
%                       w = 'Before': sum before the autocorrelation or
%                           spectrum computation.
%                       w = 'After': autocorrelation or spectrum computed
%                           for each band, and summed into a "summary".
%               mirenvelope options: 'HalfwaveCenter','Diff' (toggled on by
%                   default here),'HalfwaveDiff','Center','Smooth',
%                   'Sampling'
%               mirflux options: 'Inc','Halfwave','Complex','Median'
%       mirtempo(...,'Resonance',r) specifies the resonance curve, which
%           emphasizes the periods that are more easily perceived.
%           Possible values: 'ToiviainenSnyder' (default), 0 (toggled off)
%   Optional arguments used for the peak picking (cf. help mirpeaks)
%       mirtempo(...,'Contrast',thr): a threshold value. A given local
%           maximum will be considered as a peak if its distance with the
%           previous and successive local minima (if any) is higher than 
%           this threshold. This distance is expressed with respect to the
%           total amplitude of the autocorrelation function.
%               if no value for thr is given, the value thr=0.1 is chosen
%                   by default.
%       mirtempo(...,'Track',tr): tracks peaks along time in order to 
%           obtain a stabilized tempo curve and to limit therefore switches
%           between alternative pulsations
%               if no value for thr is given, the value tr=0.1 is chosen
%                   by default.
%
%   [t,p] = mirtempo(...) also displays the result of the signal analysis
%       leading to the tempo estimation, and shows in particular the
%       peaks corresponding to the tempo values.
            
    
        sum.key = 'Sum';
        sum.type = 'String';
        sum.choice = {'Before','After','Adjacent',0};
        sum.default = 'Before';
    option.sum = sum;
        
%% options related to mironsets:    

        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [3 .1];
    option.frame = frame;
    
        fea.type = 'String';
        fea.choice = {'Envelope','DiffEnvelope','SpectralFlux',...
                      'Pitch','Novelty'};
        fea.default = 'Envelope';
    option.fea = fea;
    
    %% options related to 'Envelope':
    
            envmeth.key = 'Method';
            envmeth.type = 'String';
            envmeth.choice = {'Filter','Spectro'};
            envmeth.default = 'Filter';
        option.envmeth = envmeth;
    
        %% options related to 'Filter':

                fb.key = 'Filterbank';
                fb.type = 'Integer';
                fb.default = 10;
            option.fb = fb;

                fbtype.key = 'FilterbankType';
                fbtype.type = 'String';
                fbtype.choice = {'Gammatone','Scheirer','Klapuri'};
                fbtype.default = 'Gammatone';
            option.fbtype = fbtype;

                ftype.key = 'FilterType';
                ftype.type = 'String';
                ftype.choice = {'IIR','HalfHann'};
                ftype.default = 'IIR';
            option.ftype = ftype;

        %% options related to 'Spectro':
        
                band.type = 'String';
                band.choice = {'Freq','Mel','Bark','Cents'};
                band.default = 'Freq';
            option.band = band;

        
            chwr.key = 'HalfwaveCenter';
            chwr.type = 'Boolean';
            chwr.default = 0;
        option.chwr = chwr;

            diff.key = 'Diff';
            diff.type = 'Boolean';
            diff.default = 1; % Different default for mirtempo
        option.diff = diff;

            diffhwr.key = 'HalfwaveDiff';
            diffhwr.type = 'Integer';
            diffhwr.default = 0;
            diffhwr.keydefault = 1;
        option.diffhwr = diffhwr;
        
            lambda.key = 'Lambda';
            lambda.type = 'Integer';
            lambda.default = 1;
        option.lambda = lambda;

            mu.key = 'Mu'; 
            mu.type = 'Integer'; 
            mu.default = 0; 
	        option.mu = mu; 
        
            log.key = 'Log';
            log.type = 'Boolean';
            log.default = 0;
        option.log = log;

            c.key = 'Center';
            c.type = 'Boolean';
            c.default = 0;
        option.c = c;

            aver.key = 'Smooth';
            aver.type = 'Integer';
            aver.default = 0;
            aver.keydefault = 30;
        option.aver = aver;

            sampling.key = 'Sampling';
            sampling.type = 'Integer';
            sampling.default = 0;
        option.sampling = sampling;

    %% options related to 'SpectralFlux'
    
            complex.key = 'Complex';
            complex.type = 'Boolean';
            complex.default = 0;
        option.complex = complex;

            inc.key = 'Inc';
            inc.type = 'Boolean';
            inc.default = 1;
        option.inc = inc;

            median.key = 'Median';
            median.type = 'Integer';
            median.number = 2;
            median.default = [.2 1.3];
        option.median = median;

            hw.key = 'Halfwave';
            hw.type = 'Boolean';
            hw.default = 1;
        option.hw = hw;                
        
        
%% options related to mirautocor:    

        aut.key = 'Autocor';
        aut.type = 'Integer';
        aut.default = 0;
        aut.keydefault = 1;
    option.aut = aut;            
    
        nw.key = 'NormalWindow';
        nw.default = 0;
    option.nw = nw;

        enh.key = 'Enhanced';
        enh.type = 'Integers';
        enh.default = 2:10;
        enh.keydefault = 2:10;
    option.enh = enh;

        r.key = 'Resonance';
        r.type = 'String';
        r.choice = {'ToiviainenSnyder','vanNoorden',0,'off','no'};
        r.default = 'ToiviainenSnyder';
    option.r = r;
    
        phase.key = 'Phase';
        phase.type = 'Boolean';
        phase.default = 0;
    option.phase = phase;

%% options related to mirspectrum:
    
        spe.key = 'Spectrum';
        spe.type = 'Integer';
        spe.default = 0;
        spe.keydefault = 1;
    option.spe = spe;

        zp.key = 'ZeroPad';
        zp.type = 'Integer';
        zp.default = 10000;
        zp.keydefault = Inf;
    option.zp = zp;
    
        prod.key = 'Prod';
        prod.type = 'Integers';
        prod.default = 0;
        prod.keydefault = 2:6;
    option.prod = prod;

    
%% options related to the peak detection

        m.key = 'Total';
        m.type = 'Integer';
        m.default = 1;
    option.m = m;
        
        thr.key = 'Threshold';
        thr.type = 'Integer';
        thr.default = 0;
    option.thr = thr;
    
        cthr.key = 'Contrast';
        cthr.type = 'Integer';
        cthr.default = 0.1;
    option.cthr = cthr;

        mi.key = 'Min';
        mi.type = 'Integer';
        mi.default = 40;
    option.mi = mi;
        
        ma.key = 'Max';
        ma.type = 'Integer';
        ma.default = 200;
    option.ma = ma;

        track.key = 'Track';
        track.type = 'Integer';
        track.keydefault = .1;
        track.default = 0;
    option.track = track;

        mem.key = 'TrackMem';
        mem.type = 'Integer';
        mem.default = 0;
        mem.keydefault = Inf;
    option.mem = mem;

        fuse.key = 'Fuse';
        fuse.type = 'Boolean';
        fuse.default = 0;
    option.fuse = fuse;

        pref.key = 'Pref';
        pref.type = 'Integer';
        pref.number = 2;
        pref.default = [0 .2];
    option.pref = pref;
            
        perio.key = 'Periodicity';
        perio.type = 'Boolean';
        perio.default = 0;
    option.perio = perio;
    
        lart.key = 'Lartillot';
        lart.type = 'Integer';
        lart.default = 0;
        lart.keydefault = .15;
    option.lart = lart;

        lart2.type = 'Integer';
        lart2.default = .15; %.2;
    option.lart2 = lart2;

            mean.key = 'Mean';
            mean.type = 'Boolean';
            mean.default = 0;
        option.mean = mean;
    
        wrap.key = 'Wrap';
        wrap.type = 'Boolean';
        wrap.default = 0;
        wrap.when = 'After';
    option.wrap = wrap;
    
specif.option = option;

varargout = mirfunction(@mirtempo,x,varargin,nargout,specif,@init,@main);


%% INIT

function [y type] = init(x,option)
if iscell(x)
    x = x{1};
end
if isamir(x,'mirscalar')
    y = x;
    return
end
if option.perio
    option.m = 3;
    option.enh = 2:10;
end
if option.track
    option.enh = 0;
end
if option.lart
    option.m = Inf;
    option.enh = 0;
    option.r = 0;
    option.mi = 24;
    option.ma = Inf;
    option.fea = 'Novelty';
    option.cthr = .05;
    %option.thr = .05;
end
if not(isamir(x,'mirautocor')) && not(isamir(x,'mirspectrum'))
    if isframed(x) && strcmpi(option.fea,'Envelope') && not(isamir(x,'mirscalar'))
        warning('WARNING IN MIRTEMPO: The input should not be already decomposed into frames.');
        disp('Suggestion: Use the ''Frame'' option instead.')
    end
    if strcmpi(option.sum,'Before')
        optionsum = 1;
    elseif strcmpi(option.sum,'Adjacent')
        optionsum = 5;
    else
        optionsum = 0;
    end
    if option.frame.length.val
        x = mironsets(x,option.fea,'Filterbank',option.fb,...
                    'FilterbankType',option.fbtype,...
                    'FilterType',option.ftype,...
                    'Sum',optionsum,'Method',option.envmeth,...
                    option.band,'Center',option.c,...
                    'HalfwaveCenter',option.chwr,'Diff',option.diff,...
                    'HalfwaveDiff',option.diffhwr,'Lambda',option.lambda,...
                    'Smooth',option.aver,'Sampling',option.sampling,...
                    'Complex',option.complex,'Inc',option.inc,...
                    'Median',option.median(1),option.median(2),...
                    'Halfwave',option.hw,'Detect',0,...
                    'Mu',option.mu,'Log',option.log,...
                    'Frame',option.frame.length.val,...
                            option.frame.length.unit,...
                            option.frame.hop.val,...
                            option.frame.hop.unit);
    else
        x = mironsets(x,option.fea,'Filterbank',option.fb,...
                    'FilterbankType',option.fbtype,...
                    'FilterType',option.ftype,...
                    'Sum',optionsum,'Method',option.envmeth,...
                    option.band,'Center',option.c,...
                    'HalfwaveCenter',option.chwr,'Diff',option.diff,...
                    'HalfwaveDiff',option.diffhwr,'Lambda',option.lambda,...
                    'Smooth',option.aver,'Sampling',option.sampling,...
                    'Complex',option.complex,'Inc',option.inc,...
                    'Median',option.median(1),option.median(2),...
                    'Halfwave',option.hw,'Detect',0,...
                    'Mu',option.mu,'Log',option.log);
    end
end
if option.aut == 0 && option.spe == 0
    option.aut = 1;
end
if isamir(x,'mirautocor') || (option.aut && not(option.spe))
    y = mirautocor(x,'Min',60/option.ma,'Max',60/option.mi,...
          'Enhanced',option.enh,...'NormalInput','coeff',...
          'Resonance',option.r,'NormalWindow',option.nw,...
          'Phase',option.phase);
elseif isamir(x,'mirspectrum') || (option.spe && not(option.aut))
    y = mirspectrum(x,'Min',option.mi/60,'Max',option.ma/60,...
                       'Prod',option.prod,...'NormalInput',...
                       'ZeroPad',option.zp,'Resonance',option.r);
elseif option.spe && option.aut
    ac = mirautocor(x,'Min',60/option.ma,'Max',60/option.mi,...
          'Enhanced',option.enh,...'NormalInput','coeff',...
          'Resonance',option.r);
    sp = mirspectrum(x,'Min',option.mi/60,'Max',option.ma/60,...
                       'Prod',option.prod,...'NormalInput',...
                       'ZeroPad',option.zp,'Resonance',option.r);
    y = ac*sp;
end
if ischar(option.sum)
    y = mirsum(y);
end
y = mirpeaks(y,'Total',option.m,'Track',option.track,...
               'TrackMem',option.mem,'Fuse',option.fuse,...
               'Pref',option.pref(1),option.pref(2),...
               'Threshold',option.thr,'Contrast',option.cthr,...
               'NoBegin','NoEnd',...
               'Normalize','Local','Order','Amplitude');
if option.phase
    y = mirautocor(y,'Phase');
end
type = {'mirscalar',mirtype(y)};            


%% MAIN

function o = main(p,option,postoption)
if iscell(p)
    p = p{1};
end
if isamir(p,'mirscalar')
    t = modif(p,postoption);
    o = {t};
    return
end
pt = get(p,'TrackPrecisePos');
track = 1;
if isempty(pt) || isempty(pt{1})
    pt = get(p,'PeakPrecisePos');
    track = 0;
end
bpm = cell(1,length(pt));
if 0 %option.lart
    meanbpm = cell(1,length(pt));
    d = get(p,'Data');
    pp = get(p,'Pos');
    ppp = get(p,'PeakPos');
    for j = 1:length(pt)
        bpm{j} = cell(1,length(pt{j}));
        meanbpm{j} = cell(1,length(pt{j}));
        for k = 1:length(pt{j})
            ptk = pt{j}{k};
            %tmpk = cell(1,size(ptk,2),size(ptk,3));
            for h = 1:size(ptk,3)
                tmpseg = [];
                %figure,hold on
                for l = 1:size(ptk,2)
                    res = 0;
                    if l == 1 || isempty(ptl)
                        % New init phase: strongest periodicity selected
                        ptl = getbpm(p,ptk{1,l,h});
                        for i = length(tmpseg):-1:1
                            for i2 = 1:size(tmpseg(i).bpms,1)
                                if abs(log(tmpseg(i).bpms(i2,end)) ...
                                        - log(ptl(1))) < option.lart
                                    res = 1;
                                    buf = tmpseg(i);
                                    tmpseg(i) = [];
                                    tmpseg(end+1) = buf;
                                    tmpseg(end).timindx(end+1) = l;
                                    tmpseg(end).bpms(i,end+1) = ptl(l);
                                    tmpseg(end).scores(i,end+1) = ...
                                        d{j}{k}(ppp{j}{k}{1,l,h}(1),l,h);
                                    break
                                end
                            end
                        end
                        if ~res
                            tmpseg(end+1).timindx = l;
                            tmpseg(end).bpms = ptl(1);
                            tmpseg(end).scores = ...
                                d{j}{k}(ppp{j}{k}{1,l,h}(1),l,h);
                        end
                        change = 0;
                        
                    else
                        ptl = getbpm(p,ptk{1,l,h});
                        for i = 1:size(tmpseg(end).bpms,1)
                            if abs(log(ptl(1)) - ...
                                    log(tmpseg(end).bpms(i,end))) ...
                                    < option.lart
                                % Continuation of strongest track
                                res = 1;
                                tmpseg(end).timindx(end+1) = l;
                                tmpseg(end).bpms(i,end+1) = ptl(1);
                                tmpseg(end).scores(i,end+1) = ...
                                    d{j}{k}(ppp{j}{k}{1,l,h}(1),l,h);
                                break
                            end
                        end
                        
                        if ~res
                            div = ptl(1)./tmpseg(end).bpms(:,end);
                            r1 = mod(div,1);
                            r2 = mod(1./div,1);
                            if ~isempty(find(r1 < option.lart2)) || ...
                                   ~isempty(find(r2 < option.lart2)) || ...
                                   ~isempty(find(r1 > 1 - option.lart2)) || ...
                                   ~isempty(find(r2 > 1 - option.lart2))
                                % Subdivision/multiple of pulse detected
                                res = 1;
                                tmpseg = newtrack(tmpseg(end),l,h,p,ptk,ptl(1),d{j}{k},pp{j}{k},ppp{j}{k},option.lart);                                
                            end
                        end
                                                
                        if ~res
                            % New track started
                            tmpo = ptl(1);
                            res = 0;
                            for i = length(tmpseg):-1:1
                                for i2 = 1:size(tmpseg(i).bpms,1)
                                    if abs(log(tmpseg(i).bpms(i2,end)) ...
                                            - log(ptl(1))) < option.lart
                                        res = 1;
                                        buf = tmpseg(i);
                                        tmpseg(i) = [];
                                        tmpseg(end+1) = buf;
                                        tmpseg(end).timindx(end+1) = l;
                                        tmpseg(end).bpms(i2,end+1) = tmpo;
                                        tmpseg(end).scores(i2,end+1) = ...
                                            d{j}{k}(ppp{j}{k}{1,l,h}(1),l,h);
                                        break
                                    end
                                end
                                if res
                                    break
                                end
                            end
                            if ~res
                                for i = length(tmpseg):-1:1
                                    for i2 = 1:size(tmpseg(i).bpms,1)
                                        div = ptl(1)./tmpseg(i).bpms(i2,end);
                                        r1 = mod(div,1);
                                        r2 = mod(1./div,1);
                                        if r1 < option.lart2 || ...
                                                r1 > 1-option.lart2 || ...
                                                r2 < option.lart2 || ...
                                                r2 > 1-option.lart2
                                            res = 1;
                                            buf = tmpseg(i);
                                            tmpseg(i) = [];
                                            tmpseg(end+1) = buf;
                                            
                                            tmpseg = newtrack(tmpseg(end),l,h,p,ptk,tmpo,d{j}{k},pp{j}{k},ppp{j}{k},option.lart);
                                            break
                                        end
                                    end
                                    if res
                                        break
                                    end
                                end
                            end
                            if ~res
                                tmpseg(end+1).timindx = l;
                                tmpseg(end).bpms = tmpo;
                                tmpseg(end).scores = ...
                                    d{j}{k}(ppp{j}{k}{1,l,h}(1),l,h);
                                %tmpseg = integrate(tmpseg,1);
                            end
                        end
                        
                        last = size(tmpseg(end).bpms,2);
                        i = 1;
                        while i <= size(tmpseg(end).bpms,1)
                            if ~tmpseg(end).bpms(i,end)
                                [tmpseg(end) res] = filltmpseg(...
                                    tmpseg(end),i,last-1,last,...
                                    ptl,d{j}{k}(:,l,h),p,...
                                    pp{j}{k}(:,l,h),ppp{j}{k}{1,l,h},...
                                    option.lart);
                                if ~res
                                    candid = tmpseg(end).scores(:,end);
                                    candid(~tmpseg(end).scores(:,end-1)) = 0;
                                    [candid best] = max(candid);
                                    if candid
                                        tmpseg(end).bpms(i,end) = ...
                                            tmpseg(end).bpms(best,end)/...
                                            tmpseg(end).bpms(best,end-1)*...
                                            tmpseg(end).bpms(i,end-1);
                                    else
                                        tmpseg(end).bpms(i,end) = ...
                                            tmpseg(end).bpms(i,end-1);
                                    end
                                end
                            end
                            
                            res = 0;
                            [tmpseg(end) res] = fusetracks(tmpseg(end),i,...
                                    size(tmpseg(end).bpms,2),option.lart,0);
                            if ~res
                                i = i+1;
                            end
                        end
                        
                    end
                    
                    %plot(tmpseg(end).timindx(end),...
                    %     60./tmpseg(end).bpms(:,end)','+-')
                    %drawnow
                end

                figure
                plot(tmpseg.timindx,60./tmpseg.bpms','+-')
                
                longest = 1;
                for i = 1:length(tmpseg)
                    if length(tmpseg(i).timindx) > ...
                            length(tmpseg(longest).timindx)
                        longest = i;
                    end
                end
                tmpseg = tmpseg(longest);
                
                [unused best] = max(mean(tmpseg.scores,2));
                bpmk = zeros(1,size(ptk,2));
                bpmk(tmpseg.timindx) = tmpseg.bpms(best,:);
                for i = 2:length(bpmk)
                    if ~bpmk(i)
                        bpmk(i) = bpmk(i-1);
                    end
                end
                bpm{j}{k} = bpmk;
            end
        end 
    end
elseif option.lart
    bpm = cell(1,length(pt));
    d = get(p,'Data');
    pp = get(p,'Pos');
    ppp = get(p,'PeakPos');
    for j = 1:length(pt)
        bpm{j} = cell(1,length(pt{j}));
        meanbpm{j} = cell(1,length(pt{j}));
        for k = 1:length(pt{j})
            ptk = pt{j}{k};
            for h = 1:size(ptk,3)
                meters = {};
                oldmeters = {};
                bpmk = zeros(1,size(ptk,2));
                currentbpmk = [];
                %figure, hold on
                for l = 1:size(ptk,2)
                    ptl = getbpm(p,ptk{1,l,h});
                    %comet = zeros(1,length(meters));
                    for i = 1:length(ptl)
                        res = 0;
                        for i2 = 1:length(meters)
                            res3 = 0;
                            bpms = [meters{i2}.lastbpm];
                            if 1
                                dist = abs(60/ptl(i) - 60./bpms);
                            else
                                norb = bpms ./ [meters{i2}.lvl];
                                nbpms = mean(norb) * [meters{i2}.lvl];
                                dist = min(abs(60/ptl(i) - 60./bpms),...
                                           abs(60/ptl(i) - 60./nbpms));
                            end
                            if 0 %comet(i2)
                                thr = option.lart * 2;
                            else
                                thr = option.lart;
                            end
                            
                            i3 = 1;
                            while i3 <= length(dist)
                                if dist(i3) < thr
                                    % Continuing an existing metrical
                                    % level.
                                    if meters{i2}(i3).timidx(end) == l
                                        % Already continued.
                                        if dist(i3) < min(abs(60./meters{i2}(i3).bpms(end) ...
                                                              - 60./bpms(i3)),...
                                                          abs(60./meters{i2}(i3).bpms(end) ...
                                                              - 60./nbpms(i3)))
                                            % New candidate is better.
                                            meters{i2}(i3).bpms(end) = ptl(i);
                                            meters{i2}(i3).score(end) = ...
                                                d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                        end
                                        res3 = [i2 i3];
                                    else
                                        if bpms(i3) > ptl(i)
                                            div = bpms(i3) / ptl(i);
                                        else
                                            div = ptl(i) / bpms(i3);
                                        end
                                        if div < 1+option.lart2
                                            if 1 % meters{i2}(i3).timidx(end) < l
                                                if res3
                                                    % Level identified to
                                                    % one already detected
                                                    meters{i2}(i3) = [];
                                                    bpms(i3) = [];
                                                    dist(i3) = [];
                                                    if ~isempty(currentbpmk) && ...
                                                            currentbpmk(1) == i2
                                                        if currentbpmk(2) == i3
                                                            currentbpmk = res3;
                                                        elseif currentbpmk(2) > i3
                                                            currentbpmk(2) = ...
                                                                currentbpmk(2) - 1;
                                                        end
                                                    end
                                                    i3 = i3 - 1;
                                                else
                                                    meters{i2}(i3).timidx(end+1) = l;
                                                    meters{i2}(i3).bpms(end+1) = ptl(i);
                                                    meters{i2}(i3).lastbpm = ptl(i);
                                                    meters{i2}(i3).score(end+1) = ...
                                                        d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                                    meters{i2}(i3).main(end+1) = ...
                                                        meters{i2}(i3).main(end);
                                                    %comet(i2) = 1;
                                                end
                                            elseif length(meters{i2}(i3).bpms) > 1 ...
                                                   && abs(meters{i2}(i3).lastbpm ...
                                                          - ptl(i)) < ...
                                                      abs(diff(meters{i2}(i3)...
                                                                    .bpms(end-1:end)))
                                                meters{i2}(i3).bpms(end) = ptl(i);
                                                meters{i2}(i3).lastbpm = ptl(i);
                                                meters{i2}(i3).score(end) = ...
                                                    d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                                meters{i2}(i3).main(end+1) = ...
                                                    meters{i2}(i3).main(end);
                                                %comet(i2) = 1;
                                            end
                                            res = 1;
                                            res3 = [i2 i3];
                                        end
                                    end
                                end
                                i3 = i3 + 1;
                            end
                        end
                        if ~res
                            % Candidate level not belonging to current
                            % metrical levels
                            i2 = 1;
                            while i2 <= length(meters)
                                bpms = [meters{i2}.lastbpm];
                                res3 = 0;
                                i3 = 1;
                                while i3 <= length(bpms)
                                    if ptl(i) > bpms(i3)
                                        % Candidate faster than stored
                                        % level
                                        div = ptl(i) / bpms(i3);
                                        rdiv = round(div);
                                        if rdiv > 1 && ...
                                                ~isempty(find(~mod(rdiv,[2 3]))) && ...
                                                (mod(div,1) < option.lart2 || ...
                                                 mod(div,1) > 1-option.lart2)
                                            % Candidate level can be
                                            % integrated in this metrical
                                            % hierarchy
                                            
                                            if res3
                                                % Level identified to one
                                                % already detected
                                                
                                                %meters{i2}(i3) = [];
                                                %bpms(i3) = [];
                                                %if ~isempty(currentbpmk) && ...
                                                %        currentbpmk(1) == i2 && ...
                                                %        currentbpmk(2) >= i3
                                                %    currentbpmk = [res3(1)...
                                                %                   res3(2)-1];
                                                %end
                                                %i3 = i3 - 1;
                                                %if res(1) == i2 && res(2)>i3
                                                %    res(2) = res(2)-1;
                                                %    res3 = res;
                                                %end
                                                
                                            elseif res
                                                % Candidate level also
                                                % integrated in other
                                                % metrical hierarchy. Both
                                                % hierarchies are fused.
                                                
                                                chcur = 0;
                                                if meters{res(1)}(res(2)).lvl > ...
                                                        meters{i2}(i3).lvl * rdiv
                                                    % Other hierarchy is
                                                    % faster than current.
                                                    meter1 = meters{i2};
                                                        % Slower hierarchy,
                                                        % which is fused to
                                                        % faster one
                                                    meter2 = meters{res(1)};
                                                        % Faster hierarchy,
                                                        % onto which slower
                                                        % one is fused
                                                    lvl1 = meter1(i3).lvl * rdiv;
                                                        % Level in slower
                                                    lvl2 = meter2(res(2)).lvl;
                                                        % Level in faster
                                                    if ~isempty(currentbpmk) && ...
                                                            currentbpmk(1) == i2
                                                        currentbpmk(1) = res(1);
                                                        chcur = 1;
                                                        % Best bpm was in
                                                        % current
                                                        % hierarchy.
                                                    end
                                                else
                                                    % Other hierarchy is
                                                    % slower than current
                                                    meter1 = meters{res(1)};
                                                    meter2 = meters{i2};
                                                    lvl1 = meter1(res(2)).lvl;
                                                    lvl2 = meter2(i3).lvl * rdiv;
                                                    if ~isempty(currentbpmk) && ...
                                                            currentbpmk(1) == res(1)
                                                        chcur = 1;
                                                        % Best bpm was in
                                                        % other hierarchy.
                                                    end
                                                end
                                                div = lvl2 / lvl1;
                                                if round(div) == div
                                                    % Faster hierarchy is
                                                    % exact multiple of
                                                    % slower one.
                                                    mult1 = round(div);
                                                    mult2 = 1;
                                                else
                                                    mult1 = lvl2;
                                                    mult2 = lvl1;
                                                end
                                                if mult2 > 1
                                                    for i4 = 1:length(meter2)
                                                        meter2(i4).lvl = meter2(i4).lvl * mult2;
                                                    end
                                                end
                                                for i4 = 1:length(meter1)
                                                    % Fusing one hierarchy
                                                    % into the other..
                                                    found = find(meter1(i4).lvl * mult1 ...
                                                                 == [meter2.lvl]);
                                                    if isempty(found)
                                                        meter2(end+1) = meter1(i4);
                                                        meter2(end).lvl = meter2(end).lvl * mult1;
                                                        if chcur == 1 && ...
                                                                i4 == currentbpmk(2)
                                                            currentbpmk(2) = length(meter2);
                                                            % Best bpm now
                                                            % pointing to
                                                            % the fused
                                                            % hierarchy.
                                                        end
                                                    end
                                                end
                                                if ~isempty(currentbpmk)
                                                    if currentbpmk(1) > i2
                                                        currentbpmk(1) = currentbpmk(1) - 1;
                                                    elseif currentbpmk(1) == i2
                                                        currentbpmk(1) = res(1);
                                                    end
                                                end
                                                meters{res(1)} = meter2;
                                                meters(i2) = [];
                                                i2 = i2 - 1;
                                                break
                                            else
                                                lvl = meters{i2}(i3).lvl ...
                                                    * round(div);
                                                l0 = find(lvl == ...
                                                          [meters{i2}.lvl]);
                                                score = d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                                if isempty(l0)
                                                    % New metrical level
                                                    meters{i2}(end+1).lvl = lvl;
                                                    meters{i2}(end).lastbpm = ptl(i);
                                                    meters{i2}(end).bpms = ptl(i);
                                                    meters{i2}(end).timidx = l;
                                                    meters{i2}(end).score = score;
                                                    %meters{i2}(end).slow.index = meters{i2}(i3).index;
                                                    %meters{i2}(end).slow.ratio = round(div);
                                                    %meters{i2}(end).fast = [];
                                                    meters{i2}(end).main = 0;
                                                    res = [i2 length(meters{i2})];
                                                    %comet(i2) = 1;
                                                elseif score > meters{i2}(l0).score
                                                    meters{i2}(l0).lastbpm = ptl(i);
                                                    meters{i2}(l0).bpms = ptl(i);
                                                    meters{i2}(l0).score = score;
                                                    res = [i2 l0];
                                                end
                                            end
                                            res3 = res;
                                        end
                                    else
                                        % Candidate slower than stored
                                        % level
                                        div = bpms(i3) / ptl(i);
                                        rdiv = round(div);
                                        if rdiv > 1 && ...
                                                ~isempty(find(~mod(rdiv,[2 3]))) && ...
                                                (mod(div,1) < option.lart2 || ...
                                                 mod(div,1) > 1-option.lart2)
                                            % Candidate level can be
                                            % integrated in this metrical
                                            % hierarchy
                                            
                                            if res3
                                                % Level identified to one
                                                % already detected
                                                
                                                %meters{i2}(i3) = [];
                                                %bpms(i3) = [];
                                                %if ~isempty(currentbpmk) && ...
                                                %        currentbpmk(1) == i2 && ...
                                                %        currentbpmk(2) >= i3
                                                %    currentbpmk = [res3(1)...
                                                %                   res3(2)-1];;
                                                %end
                                                %i3 = i3 - 1;
                                                %if res(1) == i2 && res(2)>i3
                                                %    res(2) = res(2)-1;
                                                %    res3 = res;
                                                %end
                                            else
                                                score = d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                                lvl = meters{i2}(i3).lvl ...
                                                        / round(div);
                                                if lvl == round(lvl)
                                                    l0 = find(lvl == ...
                                                              [meters{i2}.lvl]);
                                                    if isempty(l0)
                                                        % New metrical level
                                                        meters{i2}(end+1).lvl = lvl;
                                                        meters{i2}(end).lastbpm = ptl(i);
                                                        meters{i2}(end).bpms = ptl(i);
                                                        meters{i2}(end).timidx = l;
                                                        meters{i2}(end).score = score;
                                                        %meters{i2}(end).slow = [];
                                                        %meters{i2}(end).fast.index = meters{i2}(i3).index;
                                                        %meters{i2}(end).fast.ratio = round(div);
                                                        meters{i2}(end).main = 0;
                                                        res = [i2 length(meters{i2})];
                                                        %comet(i2) = 1;
                                                    elseif score > meters{i2}(l0).score
                                                        meters{i2}(l0).lastbpm = ptl(i);
                                                        meters{i2}(l0).bpms = ptl(i);
                                                        meters{i2}(l0).score = score;
                                                        res = [i2 l0];
                                                    end
                                                else
                                                    lvl = meters{i2}(i3).lvl;
                                                    for i4 = 1:length(meters{i2})
                                                        meters{i2}(i4).lvl = meters{i2}(i4).lvl * round(div);
                                                    end
                                                    meters{i2}(end+1).lvl = lvl;
                                                    meters{i2}(end).lastbpm = ptl(i);
                                                    meters{i2}(end).bpms = ptl(i);
                                                    meters{i2}(end).timidx = l;
                                                    meters{i2}(end).score = score;
                                                    meters{i2}(end).main = 0;
                                                    res = [i2 length(meters{i2})];
                                                    comet(i2) = 1;
                                                end
                                            end
                                            res3 = res;
                                        end
                                    end
                                    i3 = i3 + 1;
                                end
                                i2 = i2 + 1;
                            end
                            %for i2 = 1:length(meters)
                            %    if ~isempty(meters{i2})
                            %        [unused ord] = sort([meters{i2}.lastbpm],'descend');
                            %        if ~isempty(currentbpmk) && ...
                            %                currentbpmk(1) == i2
                            %            currentbpmk(2) = find(ord == currentbpmk(2),1);
                            %        end
                            %        meters{i2} = meters{i2}(ord);
                            %    end
                            %end
                            
                            if ~res
                                % New metrical hierarchy
                                meters{end+1}.lvl = 1;
                                meters{end}.lastbpm = ptl(i);
                                meters{end}.bpms = ptl(i);
                                meters{end}.timidx = l;
                                meters{end}.score = ...
                                    d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                meters{end}.slow = [];
                                meters{end}.fast = [];
                                meters{end}.main = 0;
                                %comet(end+1) = 1;
                            end
                        end
                    end
                    
                    for i = 2:length(meters)
                        if length(oldmeters)<i || isempty(oldmeters{i})
                            included = 1;
                            for i2 = 1:length(meters{i})
                                found = 0;
                                for i3 = 1:i-1
                                    bpms = [meters{i3}.lastbpm];
                                    norb = bpms ./ [meters{i3}.lvl];
                                    nbpms = mean(norb) * [meters{i3}.lvl];
                                    dist = min(abs(60/meters{i}(i2).lastbpm - 60./bpms),...
                                               abs(60/meters{i}(i2).lastbpm - 60./nbpms));
                                    if ~isempty(find(dist<option.lart));
                                        found = 1;
                                        break
                                    end
                                end
                                if ~found
                                    for i3 = 1:i-1
                                        for i4 = 1:length(meters{i3})
                                            ma = max(meters{i}(i2).lastbpm,...
                                                     meters{i3}(i4).lastbpm);
                                            mi = min(meters{i}(i2).lastbpm,...
                                                     meters{i3}(i4).lastbpm);
                                            div = ma / mi;
                                            rdiv = round(div);
                                            if rdiv > 1 && ...
                                                    ~isempty(find(~mod(rdiv,[2 3]))) && ...
                                                    (mod(div,1) < option.lart2 || ...
                                                     mod(div,1) > 1-option.lart2)
                                                if meters{i}(i2).lastbpm > ...
                                                        meters{i3}(i4).lastbpm
                                                    meters{i3}(end+1).lvl = ...
                                                        meters{i}(i2).lvl * rdiv;
                                                else
                                                    meters{i3}(end+1).lvl = ...
                                                        meters{i}(i2).lvl / rdiv;
                                                end
                                                meters{i3}(end).lastbpm = meters{i}(i2).lastbpm;
                                                meters{i3}(end).bpms = meters{i}(i2).bpms;
                                                meters{i3}(end).timidx = meters{i}(i2).timidx;
                                                meters{i3}(end).score = meters{i}(i2).score;
                                                meters{i3}(end).main = meters{i}(i2).main;
                                                found = 1;
                                                break
                                            end
                                        end
                                        if found
                                            break
                                        end
                                    end
                                end
                                if ~found
                                    included = 0;
                                    break
                                end
                            end
                            if included
                                % meters{i} is completely included into
                                % meters{i2}
                                meters(i) = [];
                                if length(oldmeters) >= i
                                    oldmeters(i) = [];
                                end
                                break
                            end
                        end
                    end
                    
                    mi = [];
                    for i = 1:length(meters)
                        if ~isempty(meters{i})
                            for i2 = 1:length(meters{i})
                                if meters{i}(i2).timidx(end) == l && ...
                                        (isempty(mi) || ...
                                         meters{i}(i2).score(end) > ...
                                         meters{mi(1)}(mi(2)).score(end));
                                    mi = [i i2];
                                end
                            end
                        end
                    end
                    if 0 %~isempty(mi)
                        meters{mi(1)}(mi(2)).main(end) = 1;
                        if isempty(currentbpmk) || ...
                                ~meters{currentbpmk(1)}...
                                    (currentbpmk(2)).main(end) || ...
                                meters{currentbpmk(1)}...
                                    (currentbpmk(2)).timidx(end) ~= l || ...
                                ((length(meters{mi(1)}(mi(2)).main) == 1 || ...
                                  ~meters{mi(1)}(mi(2)).main(end-1)) && ...
                                 resonance(meters{mi(1)}(mi(2)).bpms(end)) > ...
                                    resonance(meters{currentbpmk(1)}...
                                     (currentbpmk(2)).bpms(end)) * 1.2 && ...
                                 resonance(meters{mi(1)}(mi(2)).bpms(end)) * ...
                                    meters{mi(1)}(mi(2)).score(end) > ...
                                 resonance(meters{currentbpmk(1)}...
                                     (currentbpmk(2)).bpms(end)) * ...
                                    meters{currentbpmk(1)}...
                                        (currentbpmk(2)).score(end)) * 1.2
                            currentbpmk = mi;
                        end
                        bpmk(l) = meters{currentbpmk(1)}...
                                        (currentbpmk(2)).bpms(end);
                    else
                        currentbpmk = [];
                        bpmk(l) = 0;
                    end
                    
                    for i = 1:length(meters)
                        i2 = 1;
                        while i2 <= length(meters{i})
                            if meters{i}(i2).timidx(end) ~= l || ...
                                    l == size(ptk,2)
                                if length(oldmeters) < i
                                    oldmeters{i} = meters{i}(i2);
                                else
                                    oldmeters{i}(end+1) = meters{i}(i2);
                                end
                                meters{i}(i2) = [];
                                if ~isempty(currentbpmk) && ...
                                        currentbpmk(1) == i && ...
                                        i2 < currentbpmk(2)
                                    currentbpmk(2) = currentbpmk(2)-1;
                                end
                            else
                                meters{i}(i2).lastbpm = meters{i}(i2).bpms(end);
                                i2 = i2+1;
                            end
                        end
                    end
                    if ~mod(l,100)
                        l
                    end
                   %for zz = 1:length(meters{1})
                   %    plot(meters{1}(zz).timidx,60./meters{1}(zz).bpms,'+-');%,'Color',num2col(i))
                   %    drawnow
                   %end 
                end
                figure,hold on
                for i = 1:length(oldmeters)
                    irgb = shiftdim(1-num2col(i),-1);
                    mac = 0;
                    mic = 1;
                    for i2 = 1:length(oldmeters{i})
                        for i3 = 1:length(oldmeters{i}(i2).score)
                            if oldmeters{i}(i2).score(i3) > mac
                                mac = oldmeters{i}(i2).score(i3);
                            end
                            if oldmeters{i}(i2).score(i3) < mic
                                mic = oldmeters{i}(i2).score(i3);
                            end
                        end
                    end
                    micmac = mac-mic;
                    if ~micmac
                        micmac = 1;
                    end
                    for i2 = 1:length(oldmeters{i})
                        for i3 = 1:length(oldmeters{i}(i2).score)
                            rgb = ones(1,1,3);
                            rgb = ones(1,1,3) - ...
                                (oldmeters{i}(i2).score(i3) - mic) / micmac * irgb;
                            plot(oldmeters{i}(i2).timidx(i3),...
                                 60./oldmeters{i}(i2).bpms(i3),'+','Color',rgb);
                            %plot(oldmeters{i}(i2).timidx(i3),...
                            %     60./oldmeters{i}(i2).bpms(i3),rgb,...%'Color',rgb);
                            %     'MarkerSize',.5 + oldmeters{i}(i2).score(i3)*5);
                        end
                    end
                    %plot(60./bpmk,'k');
                end
            end
            
            meters = oldmeters;
            bpm{j}{k} = bpmk;
            
            if 0
                %longest = [1 1];
                meanscore = cell(1,length(oldmeters));
                maxscore = zeros(2,length(meters));
                for i = 1:length(meters)
                    meanscore{i} = zeros(1,length(meters{i}));
                    for i2 = 1:length(meters{i})
                        %if length(meters{i}(j).timidx) > ...
                        %        length(meters{longest(1)}(longest(2)).timidx)
                        %    longest = [i j];
                        %end
                        scores = zeros(size(ptk));
                        scores(meters{i}(i2).timidx) = meters{i}(i2).score;
                        meanscore{i}(i2) = mean(scores) * ...
                            resonance(mean(meters{i}(i2).bpms));
                    end
                    [maxscore(1,i) maxscore(2,i)] = max([meanscore{i}]);
                end
                %meters = meters{longest(1)}(longest(2));
                [unused best] = max(maxscore(1,:));
                bestmeter = meters{best}(maxscore(2,best));
                bpmk = zeros(1,size(ptk,2));
                bpmk(bestmeter.timidx) = bestmeter.bpms;
                for i = 2:length(bpmk)
                    if ~bpmk(i)
                        bpmk(i) = bpmk(i-1);
                    end
                end
                bpm{j}{k} = bpmk;
            end
        end
    end
else
    meanbpm = {};
    for j = 1:length(pt)
        bpm{j} = cell(1,length(pt{j}));
        for k = 1:length(pt{j})
            ptk = pt{j}{k};
            bpmk = cell(1,size(ptk,2),size(ptk,3));
            for h = 1:size(ptk,3)
                for l = 1:size(ptk,2)
                    ptl = ptk{1,l,h};
                    if isempty(ptl)
                        bpmk{1,l,h} = NaN;
                    else
                        bpmk{1,l,h} = getbpm(p,ptl);
                    end
                end
            end
            if track
                bpmk = bpmk{1};
            end
            bpm{j}{k} = bpmk;
        end 
    end
end
if option.mean
    fp = get(p,'FramePos');
    for j = 1:length(fp);
        for k = 1:length(fp{j})
            fp{j}{k} = fp{j}{k}([1 end])';
        end
    end
    t = mirscalar(p,'Data',meanbpm,'Title','Tempo','Unit','bpm','FramePos',fp);
else
    t = mirscalar(p,'Data',bpm,'Title','Tempo','Unit','bpm');
    t = modif(t,postoption);
end
o = {t,p};


function tmpseg = integrate(tmpseg,j)
res = 0;
for i = length(tmpseg)-1:-1:1
    for i2 = 1:size(tmpseg(i).bpms,1)
        if abs(log(tmpseg(i).bpms(i2,end)) - ...
               log(tmpseg(end).bpms(j,end))) < option.lart
            res = 1;
            buf = tmpseg(i);
            tmpseg(i) = [];
            for j = 1:length(buf.timindx)
                indx = find(buf.timindx(j)>tmpseg(end).timindx, 1);
                tmpseg(end).timindx(indx:end) = [buf.timindx(j) ...
                    tmpseg(end).timindx(indx:end)];
                tmpseg(end).bpms(indx:end) = [buf.bpms(j) ...
                    tmpseg(end).bpms(indx:end)];
                tmpseg(end).scores(indx:end) = [buf.scores(j) ...
                    tmpseg(end).scores(indx:end)];
            end
            break
        end
    end
    if res
        break
    end
end


function tmpseg = newtrack(tmpseg,l,h,p,ptk,tmpo,d,pp,ppp,param)
tmpseg.timindx(end+1) = l;
tmpseg.bpms(end+1,end+1) = tmpo;
tmpseg.scores(end+1,end+1) = d(ppp{1,l,h}(1),l,h);

new = size(tmpseg.bpms,1);
nl = size(tmpseg.bpms,2);                                
for i = nl-1:-1:1
    [tmpseg res] = filltmpseg(tmpseg,new,i+1,i,getbpm(p,ptk{1,i,h}),...
                               d(:,i,h),p,pp(:,i,h),ppp{1,i,h},param);
    if ~res
        break
    end
    
    %plot(tmpseg.timindx(i),60/tmpseg.bpms(new,i)','o--')
    %drawnow
    
    [tmpseg res] = fusetracks(tmpseg,new,i,param,1);
    if res
        break
    end
end


function [tmpseg res] = filltmpseg(tmpseg,i,from,to,ptl,d,p,pp,ppp,param)
res = 0;
tmpi = tmpseg.bpms(i,from);
for i2 = 2:length(ptl)
    if abs(log(ptl(i2)) - log(tmpi)) < param
        % Continuation of secondary track
        res = 1;
        tmpseg.bpms(i,to) = ptl(i2);
        tmpseg.scores(i,to) = d(ppp(i2));
        break
    end
end
if ~res
    % Track not yet detected is continued if
    % there is energy left around that frequency.
    ps = find(pp > getpos(p,tmpi),1);
    md = min(d);
    dn = (d-md)/(max(d)-md);
    if dn(ps) > .05
        res = 1;
        tmpseg.bpms(i,to) = tmpseg.bpms(i,from);
        tmpseg.scores(i,to) = d(ps);
    end
end
                                

function [tmpseg res] = fusetracks(tmpseg,i,l,param,decr)
res = 0;
return
for i2 = 1:i-1
    if abs(log(tmpseg.bpms(i,l)) - log(tmpseg.bpms(i2,l))) < param/2
        if decr
            tmpseg.bpms(i,1:l-1) = tmpseg.bpms(i2,1:l-1);
            tmpseg.scores(i,1:l-1) = tmpseg.scores(i2,1:l-1);
        end
        [unused best] = min(mean(tmpseg.scores([i i2],:),2));
        switch best
            case 1
                i0 = i2;
            case 2
                i0 = i;
        end
        tmpseg.bpms(i0,:) = [];
        tmpseg.scores(i0,:) = [];
        res = 1;
        break
    end
end


function r = resonance(bpm)
r = 1 - 0.25*(log2(max(60/bpm,1e-12)/.5)).^2;
                            
                            
function bpm = getbpm(p,ptl)
if isa(p,'mirautocor') && not(get(p,'FreqDomain'))
    bpm = 60./ptl;
else
    bpm = ptl*60;
end


function ptl = getpos(p,bpm)
if isa(p,'mirautocor') && not(get(p,'FreqDomain'))
    ptl = 60./bpm;
else
    ptl = bpm/60;
end


function t = modif(t,option)
if option.wrap
    d = get(t,'Data');
    for i = 1:length(d)
        for j = 1:length(d{i})
            dij = zeros(1,length(d{i}{j}));
            r = NaN;
            for k = 1:length(d{i}{j})
                if isnan(r) || abs(diff(log(d{i}{j}([k-1 k])))) > .2
                    r = ceil(log2(d{i}{j}(k)) - log2(60)) - 1;
                end
                if r>0
                    dij(k) = d{i}{j}(k)*2^(-r);
                else
                    dij(k) = d{i}{j}(k);
                end
            end
            d{i}{j} = dij;
        end
    end
    t = set(t,'Data',d);
end