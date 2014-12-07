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
%   mirtempo(..., ?Metre?) tracks tempo by building a hierarchical metrical
%       structure (using mirmetre). This enables to find coherent metrical 
%       levels leading to a continuous tempo curve.
%   When the ?Metre? option is used for academic research, please cite the 
%       following publication:
%   Lartillot, O., Cereghetti, D., Eliard, K., Trost, W. J., Rappaz, M.-A.,
%       Grandjean, D., "Estimating tempo and metrical features by tracking 
%       the whole metrical hierarchy", 3rd International Conference on 
%       Music & Emotion, Jyväskylä, 2013.
%
%   mirtempo(..., ?Change?) computes the difference between successive 
%       values of the tempo curve. Tempo change is expressed independently 
%       from the choice of a metrical level by computing the ratio of tempo 
%       values between successive frames, and is expressed in logarithmic 
%       scale (base 2), so that no tempo change gives a value of 0, 
%       increase of tempo gives positive value, and decrease of tempo gives
%       negative value.
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
        r.choice = {'ToiviainenSnyder','vanNoorden',0,'off','no','New'};
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
    
        metre.key = 'Metre';
        metre.type = 'Integer';
        metre.default = 0;
        metre.keydefault = 2;
    option.metre = metre;
    
        minres.key = 'MinRes';
        minres.type = 'Integer';
        minres.default = 1; %.1;
    option.minres = minres;
    
        change.key = 'Change'; %% Does not work without 'Metre' so far..
        change.type = 'Boolean';
        change.default = 0;
        change.when = 'After';
    option.change = change;
    
        fill.key = 'Fill';
        fill.type = 'Boolean';
        fill.default = 1;
    option.fill = fill;
    
        mean.key = 'Mean';
        mean.type = 'Boolean';
        mean.default = 0;
    option.mean = mean;
        
specif.option = option;

varargout = mirfunction(@mirtempo,x,varargin,nargout,specif,@init,@main);


%% INIT

function [y type] = init(x,option)
if iscell(x)
    x = x{1};
end
if isamir(x,'mirscalar') || isamir(x,'mirmetre')
    y = x;
    type = {'mirscalar',mirtype(y)};            
    return
end

if option.metre
    if ~option.frame.length.val
        option.frame.length.val = 5;
        option.frame.length.unit = 's';
        option.frame.hop.val = .05;
        option.frame.hop.unit = '/1';
    end
    y = mirmetre(x,'Frame',option.frame.length.val,...
                           option.frame.length.unit,...
                           option.frame.hop.val,...
                           option.frame.hop.unit,...
                           'MinRes',option.minres);
else

    if option.perio
        option.m = 3;
        option.enh = 2:10;
    end
    if option.track
        option.enh = 0;
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
end

type = {'mirscalar',mirtype(y)};            


%% MAIN

function o = main(p,option,postoption)
if iscell(p)
    p = p{1};
end
    
if isamir(p,'mirscalar') 
    o = modif(p,postoption);
    return
end

if isa(p,'mirmetre')
    d = get(p,'Data');
    fp = get(p,'FramePos');
    g = get(p,'Globpm');
    bpm = cell(1,length(d));
    for j = 1:length(d)
        bpm{j} = cell(1,length(d{j}));
        for k = 1:length(d{j})
            if option.fill
                bpm{j}{k} = NaN(option.m,size(fp{j}{k},2));
            end
            scork = zeros(1,size(fp{j}{k},2));
            bestk = NaN(1,size(fp{j}{k},2));
            errok = Inf(1,size(fp{j}{k},2));
            spans = zeros(1,length(d{j}{k}));
            for i = 1:length(d{j}{k})
                % For each metrical hierarchy...
                scori = zeros(length(d{j}{k}{i}),size(fp{j}{k},2));
                
                [unused ll] = sort([d{j}{k}{i}.lvl]);%,'descend');
                d{j}{k}{i} = d{j}{k}{i}(ll);
                for l = 1:length(d{j}{k}{i})
                    % For each metrical level...
                    scori(l,d{j}{k}{i}(l).timidx) = d{j}{k}{i}(l).score;...
%                        .* max(0, 1 - abs(60./d{j}{k}{i}(l).bpms ...
%                                          - 60./d{j}{k}{i}(l).globpms) * 10);
                    
                    %for l2 = 1:l-1
                    %    if 0 && ~mod(d{j}{k}{i}(l2).lvl,d{j}{k}{i}(l).lvl) ...
                    %            ...d{j}{k}{i}(l2).lvl == 2*d{j}{k}{i}(l).lvl ...
                    %            && 60/max(d{j}{k}{i}(l2).bpms)<2.3
                    %        scori(l,:) = min(scori(l,:),scori(l2,:));
                    %    end                                       
                        %if ~mod(d{j}{k}{i}(l2).lvl,d{j}{k}{i}(l).lvl)
                        %    [unused ia ib] = intersect(d{j}{k}{i}(l).timidx,...
                        %                               d{j}{k}{i}(l2).timidx);
                        %    iab = find(scori{l}(ia));
                        %    scori{l}(ia(iab)) = ...
                        %        max(scori{l}(ia(iab)), ...
                        %            d{j}{k}{i}(l2).score(ib(iab)));
                        %end
                    %end
                    
                    idx = find(scori(l,:) > scork);% | ...
                               %(scori(l,:) == scork & ...
                               % abs(d{j}{k}{i}(l).bpms - ...
                               %     d{j}{k}{i}(l).globpms) < ...
                               %         errok(d{j}{k}{i}(l).timidx)));
                    scork(idx) = scori(l,idx);
                    errok(idx) = 0; %abs(d{j}{k}{i}(l).bpms(idx) - ...
                                    %d{j}{k}{i}(l).globpms(idx));
                    bestk(idx) = i;
                    spans(i) = max(spans(i),length(idx));
                end
                
                sd = zeros(1,length(d{j}{k}{i}));
                for l = 1:length(d{j}{k}{i})
                    sd(l) = sum(scori(l,:));% * resonance(mb);
                end
                
                %figure,plot([d{j}{k}{i}.lvl],sd,'+-');
                [unused bests] = sort(sd,'descend');
                
                if 0
                    best = bests(1:min(option.m,length(bests)));

                else
                    if option.metre == 1
                        metstr = struct;
                        metstr.lvl = d{j}{k}{i}(bests(1)).lvl;
                        metstr.score = score(d{j}{k}{i},sd,bests(1));
                        metstr.indx = bests(1);
                        for h = 2:length(bests)
                            % For each periodicity, from best to worst
                            found = 0;
                            for l = 1:length(metstr)
                                if mod(max(metstr(l).lvl,...
                                           d{j}{k}{i}(bests(h)).lvl),...
                                       min(metstr(l).lvl,...
                                           d{j}{k}{i}(bests(h)).lvl))
                                       found = 1;
                                       break
                                end
                            end
                            if ~found
                                metstr(end+1).lvl = d{j}{k}{i}(bests(h)).lvl;
                                metstr(end).score = score(d{j}{k}{i},sd,bests(h));
                                metstr(end).indx = bests(h);
                            end
                        end
                        bests = [metstr.indx];
                    else
                        metstr = {};
                        proc = [];
                        for h = 1:length(bests)
                            % For each periodicity, from best to worst
                            if ismember(bests(h),proc)
                                continue
                            end

                            metstr2 = struct;
                            metstr2.lvl = d{j}{k}{i}(bests(h)).lvl;
                            metstr2.score = score(d{j}{k}{i},sd,bests(h));
                            metstr2.indx = bests(h);

                            %for l = 1:bests(h)-1
                            %    if ~mod(d{j}{k}{i}(bests(h)).lvl,...
                            %            d{j}{k}{i}(l).lvl)
                            %        % Subharmonic, integrated in the hierarchy
                            %        %proc(end+1) = l;
                            %        metstr2(end+1).lvl = d{j}{k}{i}(l).lvl;
                            %        metstr2(end).score = score(d{j}{k}{i},sd,l);
                            %        metstr2(end).indx = l;
                            %    end
                            %end

                            metstr2 = recurs1(metstr2,d{j}{k}{i},sd,bests(h));

                            metstr = [metstr recurs2(metstr2,d{j}{k}{i},sd,bests(h))];
                        end
                        sumscore = zeros(1,length(metstr));
                        for h = 1:length(metstr)
                            sumscore(h) = sum([metstr{h}.score]);
                        end
                        [unused best] = max(sumscore);
                        bests = [metstr{best}.indx];
                    end
                    
                    sd1 = zeros(1,length(bests));
                    for l = 1:length(d{j}{k}{i})
                        d{j}{k}{i}(l).function = [];
                    end
                    for l = 1:length(bests)
                        d{j}{k}{i}(bests(l)).function = [1 1];
                        mb = median(d{j}{k}{i}(bests(l)).bpms);
                        sd1(l) = sd(bests(l));
                        if ischar(option.r)
                            sd1(l) = sd1(l) * resonance(mb,option.r);
                        end
                    end
                    [unused best] = sort(sd1,'descend');
                    best = bests(best(1:min(option.m,length(best))));
                
                end
                
                if 0
                    %d{j}{k}{i}(bests(1)).lvl;
                    stru = bests(1);
                    for l = 1:length(d{j}{k}{i})
                        if l ~= bests(1) && ...
                                (~mod(d{j}{k}{i}(bests(1)).lvl,d{j}{k}{i}(l).lvl) || ...
                                 ~mod(d{j}{k}{i}(l).lvl,d{j}{k}{i}(bests(1)).lvl))
                             stru(end+1) = l;
                        end
                    end
                    sd1 = zeros(1,length(stru));
                    for l = 1:length(stru)
                        mb = median(d{j}{k}{i}(stru(l)).bpms);
                        sd1(l) = sd(stru(l)) * resonance(mb);
                    end
                    [unused best1] = sort(sd1,'descend');
                    best = stru(best1(1));

                    if option.m == 2 && length(bests) > 1
                        if l == 1
                            stru = bests(2);
                            for l = 1:length(d{j}{k}{i})
                                if l ~= bests(1) && l ~= bests(2) && ...
                                        (~mod(d{j}{k}{i}(bests(2)).lvl,d{j}{k}{i}(l).lvl) || ...
                                         ~mod(d{j}{k}{i}(l).lvl,d{j}{k}{i}(bests(2)).lvl))
                                     stru(end+1) = l;
                                end
                            end
                            sd2 = zeros(1,length(stru));
                            for l = 1:length(stru)
                                mb = median(d{j}{k}{i}(stru(l)).bpms);
                                sd2(l) = sd(stru(l)) * resonance(mb);
                            end
                            [unused best2] = max(sd2);
                            best(2) = stru(best2);
                        else
                            best(2) = stru(best1(2));
                        end
                    end
                end
                
                idx = find(bestk == i);
                if isempty(idx)
                    continue
                end
                % Detecting contiguous region where one single hierarchy is
                % dominant
                idx0 = idx(find(diff([-Inf idx]) > 1)); % Starting times
                idx1 = [idx0(2:end)-1 idx(end)]; % Ending times
                
                for l = 1:length(idx0)
                    if idx1(l) - idx0(l) < 1
                        continue
                    end
                    %if 1 || idx0(l) == 1 || isnan(bestk(idx0(l)-1))
                        
                    if option.fill || spans(i) == max(spans(1:i))
                        
                        if option.fill
                            nl = option.m;
                        else
                            nl = min(option.m,length(best));
                            bpm{j}{k} = NaN(nl,size(fp{j}{k},2));
                        end
                        idx1l = min(idx1(l),size(g{j}{k}(i,:),2));
                        for h = 1:nl
                            bpm{j}{k}(h,idx0(l):idx1l) = ...
                                g{j}{k}(i,idx0(l):idx1l) ...
                                    / d{j}{k}{i}(best(h)).lvl;
                        end
                    end
                            
                    %else
                    %    dist = abs(bpm{j}{k}(idx0(l)-1) - ...
                    %               g{j}{k}(i,idx0(l)) ...
                    %                    ./ [d{j}{k}{i}.lvl]);
                    %    [unused bestl] = min(dist);
                    %    bpm{j}{k}(idx0(l):idx1(l)) = ...
                    %        g{j}{k}(i,idx0(l):idx1(l)) ...
                    %            / d{j}{k}{i}(bestl).lvl;
                    %end
                end
                
                %bpm{j}{k}(idx) = g{j}{k}(i,idx) / d{j}{k}{i}(best).lvl;
                
                bestlvl = d{j}{k}{i}(best(1)).lvl;
                for l = 1:length(d{j}{k}{i})
                    d{j}{k}{i}(l).lvl = d{j}{k}{i}(l).lvl / bestlvl;
                end
            end
        end 
    end
    t = mirscalar(p,'Data',bpm,'Title','Tempo','Unit','bpm');
    t = modif(t,postoption);
    p = set(p,'Data',d);
    o = {t,p};
    return
end

pt = get(p,'TrackPrecisePos');
track = 1;
if isempty(pt) || isempty(pt{1})
    pt = get(p,'PeakPrecisePos');
    track = 0;
end
bpm = cell(1,length(pt));
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


function s = score(d,sd,l)
if 1
    s = sd(l);
else
    s = 0;
    n = 0;
    for l2 = 1:length(d)
        if ~mod(d(l2).lvl,d(l).lvl)
            s = s + sd(l2);
            n = n + 1;
        end
    end
    s = s / n;
end


function metstr2 = recurs1(metstr1,d,sd,i)
metstr2 = {};
for j = i-1:-1:1
    if ~mod(d(i).lvl,d(j).lvl)
        metstr3 = metstr1;
        metstr3(end+1).lvl = d(j).lvl;
        metstr3(end).score = score(d,sd,j);
        metstr3(end).indx = j;
        metstr2 = [metstr2 recurs1(metstr3,d,sd,j)];
    end
end
if isempty(metstr2)
    metstr2 = {metstr1};
end


function metstr2 = recurs2(metstr1,d,sd,i)
metstr2 = {};
for j = i+1:length(d)
    if ~mod(d(j).lvl,d(i).lvl)
        metstr3 = metstr1;
        for k = 1:length(metstr3)
            metstr3{k}(end+1).lvl = d(j).lvl;
            metstr3{k}(end).score = score(d,sd,j);
            metstr3{k}(end).indx = j;
        end
        metstr2 = [metstr2 recurs2(metstr3,d,sd,j)];
    end
end
if isempty(metstr2)
    metstr2 = metstr1;
end


function r = resonance(bpm,type)
if bpm > 120
    if strcmpi(type,'New')
        bpm = (360 - bpm)/2;
    else
        bpm = 240 - bpm;
    end
end
r = max(1 - 0.25*(log2(max(60/bpm,1e-12)/.5)).^2, 0);
                            
                            
function bpm = getbpm(p,ptl)
if isa(p,'mirautocor') && not(get(p,'FreqDomain'))
    bpm = 60./ptl;
else
    bpm = ptl*60;
end


function t = modif(t,option)
if option.change
    d = get(t,'Data');
    for i = 1:length(d)
        for j = 1:length(d{i})
            if iscell(d{i}{j})
                % 'Metre' option not used here. 'Change' should not be
                % performed.
                return
            end
            d{i}{j} = [NaN diff(log2(d{i}{j}))];
        end
    end
    t = set(t,'Data',d,'Title','Tempo Change','Unit','log2(bpm)');
end