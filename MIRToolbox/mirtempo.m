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
        
        thr.key = 'Contrast';
        thr.type = 'Integer';
        thr.default = 0.1;
    option.thr = thr;

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
        lart.keydefault = .2;
    option.lart = lart;
    
    
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
               'Contrast',option.thr,'NoBegin','NoEnd',...
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
if option.lart
    d = get(p,'Data');
    pp = get(p,'Pos');
    for j = 1:length(pt)
        bpm{j} = cell(1,length(pt{j}));
        for k = 1:length(pt{j})
            ptk = pt{j}{k};
            tmpk = cell(1,size(ptk,2),size(ptk,3));
            bpm{j}{k} = zeros(1,size(ptk,2),size(ptk,3));
            changes = [];
            stdl = [];
            trsl = [];
            for h = 1:size(ptk,3)
                for l = 1:size(ptk,2)
                    ptl = getbpm(p,ptk{1,l,h});
                    
                    if isempty(ptl)
                        % No peak detected
                        
                        if ~isempty(stdl)
                            % Prolongation of previous peak
                            tmpo = bpm{j}{k}(1,stdl,h);
                            tmpk{1,l,h} = tmpk{1,stdl,h};
                            
                        else
                            % or frame leaved empty
                            tmpo = NaN;
                        end
                        
                        changes(end+1) = l;
                        
                        if ~isempty(trsl) && l-trsl<5
                            % Hypothetical phase aborted
                            bpm{j}{k}(1,trsl:l-1,h) = bpm{j}{k}(1,trsl-1,h);
                            for i = trsl:l-1
                                tmpk{1,i,h} = tmpk{1,trsl-1,h};
                            end
                        end
                        
                    elseif l == 1 || isempty(tmpk{1,l-1,h})
                        % New init phase: highest peak selected
                        tmpk{1,l,h} = 1;
                        tmpo = ptl(1);
                        changes(end+1) = l;
                        
                    else
                        res = 0;

                        for i = 1:length(tmpk{1,l-1,h})
                            if abs(log(ptl(1)*tmpk{1,l-1,h}(i)) - log(tmpo)) < option.lart
                                % Continuation of current track
                                tmpk{1,l,h}(i) = round(tmpo/ptl(1));
                                res = 1;
                                tmpo = ptl(1)*tmpk{1,l,h}(i);
                                break
                            end
                        end
                        
                        if ~res
                            if ptl(1)>tmpo
                                div = ptl(1)/tmpo;
                                r = mod(div,1);
                                if r<option.lart || r>1-option.lart
                                    % Subdivision of pulse detected
                                    % Track transformed and continued
                                    for g = changes(end):l-1
                                        tmpk{1,g,h} = tmpk{1,g,h}...
                                                        * round(div);
                                    end
                                    tmpk{1,l,h}(i+1) = 1;%%%
                                    bpm{j}{k}(1,changes(end):l-1,h) = ...
                                        bpm{j}{k}(1,changes(end):l-1,h)...
                                            * round(div);
                                    res = 1;
                                    tmpo = ptl(1);
                                end
                                
                            else
                                div = tmpo/ptl(1);
                                if round(div)>1
                                    % Multiple of pulse detected
                                    % Track continued
                                    r = mod(div,1);
                                    if r<option.lart || r>1-option.lart
                                        tmpk{1,l,h}(i+1) = round(div);
                                        res = 1;
                                        tmpo = ptl(1)*tmpk{1,l,h}(i+1);
                                    end
                                end
                            end
                        end
                        
                        if ~res && ~isempty(trsl)
                            % Same 3 steps but tried on hypothetical track
                            tmpol = bpm{j}{k}(1,trsl-1,h);
                            for i = 1:length(tmpk{1,trsl-1,h})
                                if abs(log(ptl(1)*tmpk{1,trsl-1,h}(i)) ...
                                        - log(tmpol)) ...
                                            < option.lart
                                    % Continuation of hypothetical track
                                    tmpk{1,l,h}(i) = round(tmpol/ptl(1));
                                    res = 1;
                                    tmpo = ptl(1)*tmpk{1,l,h}(i);
                                    break
                                end
                            end
                            
                            if ~res
                                if ptl(1)>tmpol
                                    div = ptl(1)/tmpol;
                                    r = mod(div,1);
                                    if r<option.lart || r>1-option.lart
                                        % Subdivision of pulse detected
                                        % Hypothetical track transformed and continued
                                        for g = changes(end):l-1
                                            tmpk{1,g,h} = tmpk{1,g,h}...
                                                            * round(div);
                                        end
                                        tmpk{1,l,h}(i+1) = 1;
                                        bpm{j}{k}(1,changes(end):trsl-1,h) = ...
                                            bpm{j}{k}(1,changes(end):trsl-1,h)...
                                                * round(div);
                                        res = 1;
                                        tmpo = ptl(1);
                                    end
                                    
                                else
                                    div = tmpol/ptl(1);
                                    r = mod(div,1);
                                    if r<option.lart || r>1-option.lart
                                        % Multiple of pulse detected
                                        % Hypothetical track continued
                                        tmpk{1,l,h}(i+1) = round(div);
                                        res = 1;
                                        tmpo = ptl(1)*tmpk{1,l,h}(i+1);
                                    end
                                end
                            end
                            
                            if res
                                % Hypothetical track turned into confirmed
                                % track
                                bpm{j}{k}(1,trsl:l-1,h) = ...
                                    bpm{j}{k}(1,trsl-1,h);
                                for g = trsl:l-1
                                    tmpk{1,g,h} = tmpk{1,trsl-1,h};
                                end
                                trsl = [];
                                stdl = l;
                            end
                        end
                        
                        if res
                            if l-changes(end)>5
                                % Hypothetical track aborted
                                stdl = l;
                                trsl = [];
                            end
                            
                        else
                            if ~isempty(stdl) && l-stdl > 10
                                % After long stop, previous peak forgotten
                                stdl = [];
                            end                            
                            
                            tmpo = ptl(1);
                            if ptl(1)<50
                                % Very slow pulses subdivided
                                for i = 2:length(ptl)
                                    if ptl(i)<180
                                        div = ptl(i)/ptl(1);
                                        if round(div)>1
                                            r = mod(div,1);
                                            if r<option.lart || ...
                                                    r>1-option.lart
                                                tmpo = ptl(i);
                                                break
                                            end
                                        end
                                    end
                                end
                            end
                            
                            tmpk{1,l,h} = 1;
                            changes(end+1) = l;
                            
                            if ~isempty(trsl) && l-trsl<5
                                % Hypothetical phase aborted
                                bpm{j}{k}(1,trsl:l-1,h) = bpm{j}{k}(1,trsl-1,h);
                                for i = trsl:l-1
                                    tmpk{1,i,h} = tmpk{1,trsl-1,h};
                                end
                            end
                            trsl = l;
                        end
                        
                        i = 1;
                        while i <= length(tmpk{1,l,h})
                            if ~tmpk{1,l,h}(i)
                                % Previous track not detected is continued if
                                % there is energy left around that frequency.
                                tmpi = tmpo/tmpk{1,l-1,h}(i);
                                ps = find(pp{j}{k}(:,l,h)...
                                            > getpos(p,tmpi),1);
                                if d{j}{k}(ps,l,h) > .1
                                    tmpk{1,l,h}(i) = tmpk{1,l-1,h}(i);
                                else
                                    tmpk{1,l,h}(i) = [];
                                    i = i-1;
                                end
                            end
                            i = i+1;
                        end
                    end
                    bpm{j}{k}(1,l,h) = tmpo;
                end
                
                r = Inf;
                newbpm = zeros(1,size(bpm{j}{k},2));
                for l = 1:size(bpm{j}{k},2)
                    if isinf(r) || abs(diff(log(bpm{j}{k}([l-1 l])))) > .2
                        if (bpm{j}{k}(l) < 50 || bpm{j}{k}(l) > 150)
                            ptl = getbpm(p,ptk{1,l,h});
                            found = 0;
                            for i = 1:length(ptl)
                                if (ptl(i) > 50 && ptl(i) < 160)
                                    if bpm{j}{k}(l) > ptl(i)
                                        div = bpm{j}{k}(l) / ptl(i);
                                        rm = mod(div,1);
                                        if 1 %rm < option.lart || ...
                                                %rm > 1-option.lart
                                            r = div;
                                            found = 1;
                                            break
                                        end
                                    else
                                        div = ptl(i) / bpm{j}{k}(l);
                                        rm = mod(div,1);
                                        if 1 %rm < option.lart || ...
                                                %rm > 1-option.lart
                                            r = 1 / div;
                                            found = 1;
                                            break
                                        end
                                    end
                                end
                            end
                            if ~found
                                r = Inf;
                            end
                        else
                            r = 1;
                        end
                    end
                    newbpm(l) = bpm{j}{k}(l) / r;
                end
                bpm{j}{k} = newbpm;
            end
        end 
    end
else
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
t = mirscalar(p,'Data',bpm,'Title','Tempo','Unit','bpm');
t = modif(t,postoption);
o = {t,p};


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