function varargout = mironsets(x,varargin)
%   o = mironsets(x) shows a temporal curve where peaks relate to the
%       temporal position of events, and estimates those event time 
%       positions.
%   Optional arguments:
%       mironsets(...,f) selects the strategy for the computation of the
%           onset detection function.
%           f = 'Envelope': Envelope of the audio signal. (Default choice).
%           With two methods for envelope extraction:
%               mironsets(...,'Spectro') (Default):
%                   mironsets(...,'SpectroFrame',fl,fh) species the frame
%                       length fl (in s.) and the hop factor fh (as a value
%                       between 0 and 1)
%                       Default values: fl = .1 s., fh = .1
%                    the frequency reassigment method can be specified:
%                    'Freq' (default), 'Mel', 'Bark' or 'Cents' (cf. mirspectrum).
%               mironsets(...,'Filter'):
%                   mironsets(...,'Filterbank',nc) specifies a preliminary
%                       filterbank decomposition into nc channels. If nc = 0,
%                       no decomposition is performed.
%                       Default value: 40.
%                   mironsets(...,'FilterbankType',ft) specifies the type of
%                       filterbank (see mirfilterbank).
%                       Default value: 'Gammatone';
%                   Options associated to the mirenvelope function can be
%                       passed here as well (see help mirenvelope):
%                      'FilterType','Tau','PreDecim'
%               mironsets(...,'Sum','no') does not sum back the channels at
%                   the end of the computation. The resulting onset curve
%                   remains therefore decomposed into several channels.
%               Options associated to the mirenvelope function can be
%                   passed here as well (see help mirenvelope):
%                   'HalfwaveCenter','Diff','HalfwaveDiff','Center',
%                   'Smooth', 'Sampling','Log','Power','Lambda',
%                  ,'PostDecim','UpSample'
%           f = 'SpectralFlux': Spectral flux of the audio signal.
%               Options associated to the mirflux function can be
%               passed here as well (see help mirflux):
%                   'Inc' (toggled on by default here),
%                   'Halfwave' (toggled on by default here),
%                   'Complex' (toggled off by default),
%                   'Median' (toggled on by default here)
%           f = 'Emerge': is an improved version of the 'SpectralFlux'
%               method that is able to detect more notes and in the same 
%               time ignore the spectral variation produced by vibrato.
%%%%
%   When the 'Emerge' method is used for academic research, please cite the 
%       following publication:
%   Lartillot, O., Cereghetti, D., Eliard, K., Trost, W. J., Rappaz, M.-A.,
%       Grandjean, D., "Estimating tempo and metrical features by tracking 
%       the whole metrical hierarchy", 3rd International Conference on 
%       Music & Emotion, Jyv?skyl?, 2013.
%%%%
%           f = 'Pitch ':computes a frame-decomposed autocorrelation function ,
%                of same default characteristics than those returned
%                by mirpitch, with however a range of frequencies set by 
%                the following options:
%                   'Min' (set by default to 30 Hz),
%                   'Max' (set by default to 1000 Hz),
%                and subsequently computes the novelty curve of the 
%                resulting similatrix matrix.
%               Option associated to the mirnovelty function can be
%               passed here as well (see help mirnovelty):
%                   'KernelSize' (set by default to 32 samples)
%       mironsets(...,'Detect',d) toggles on or off the onset detection, 
%           which is based on the onset detection function.
%           (By default toggled on.)
%           Option associated to the mirpeaks function can be specified as
%               well:
%               'Contrast' with default value c = .01
%               'Threshold' with default value t = 0
%               'Single' detects only the highest peak.
%       mironsets(...,'Attack') (or 'Attacks') detects attack phases.
%       mironsets(...,'Decay') (or 'Decays') detects decay phases.
%       mironsets(...,'Frame',...) decomposes into frames, with default frame
%           length 3 seconds and hop factor .1
%   Preselected onset detection models:
%       mironsets(...,'Scheirer') corresponds to (Scheirer, 1998):
%           mironsets(...,'FilterBankType','Scheirer',...
%                         'FilterType','HalfHann','Sampling',200,...
%                         'HalfWaveDiff','Sum',0,'Detect',0)
%       mironsets(...,'Klapuri99') corresponds to most of (Klapuri, 1999).
        
%% options related to 'Envelope':

        env.key = 'Envelope';
        env.type = 'Boolean';
        env.default = NaN;
    option.env = env;

        envmethod.key = 'Method'; % optional
        envmethod.type = 'Boolean';
    option.envmethod = envmethod;
    
        envmeth.type = 'String';
        envmeth.choice = {'Filter','Spectro'};
        envmeth.default = 'Spectro';
    option.envmeth = envmeth;
 
%%      options related to 'Filter':

            filter.key = 'FilterType';
            filter.type = 'String';
            filter.choice = {'IIR','HalfHann','Butter'};
            filter.default = NaN;
        option.filter = filter;

            tau.key = 'Tau';
            tau.type = 'Integer';
            tau.default = .02;
        option.tau = tau;

            cutoff.key = 'CutOff';
            cutoff.type = 'Integer';
            cutoff.default = 37;
    option.cutoff = cutoff;
    
            fb.key = {'Filterbank','NbChannels'};
            fb.type = 'Integer';
            fb.default = NaN;
        option.fb = fb;

            filtertype.key = 'FilterbankType';
            filtertype.type = 'String';
            %filtertype.choice = {'Gammatone','2Channels','Scheirer','Klapuri'};
            filtertype.default = 'Gammatone';
        option.filtertype = filtertype;

            decim.key = {'Decim','PreDecim'};
            decim.type = 'Integer';
            decim.default = 0;
        option.decim = decim;

            hilb.key = {'Hilbert'};
            hilb.type = 'Boolean';
            hilb.default = NaN;
        option.hilb = hilb;        
        
%%      options related to 'Spectro':

            band.type = 'String';
            band.choice = {'Freq','Mel','Bark','Cents'};
            band.default = 'Freq';
        option.band = band;
        
            specframe.key = 'SpectroFrame';
            specframe.type = 'Integer';
            specframe.number = 2;
            specframe.default = NaN;
        option.specframe = specframe;
        
            presilence.key = 'PreSilence';
            presilence.type = 'Boolean';
            presilence.default = 0;
        option.presilence = presilence;

            postsilence.key = 'PostSilence';
            postsilence.type = 'Boolean';
            postsilence.default = 0;
        option.postsilence = postsilence;
        
            powerspectrum.key = 'PowerSpectrum';
            powerspectrum.type = 'Boolean';
            powerspectrum.default = 1;
        option.powerspectrum = powerspectrum;        

            timesmooth.key = 'TimeSmooth';
            timesmooth.type = 'Boolean';
            timesmooth.default = 0;
            timesmooth.keydefault = 10;
        option.timesmooth = timesmooth;        

            terhardt.key = 'Terhardt';
            terhardt.type = 'Boolean';
            terhardt.default = 0;
        option.terhardt = terhardt;

        sum.key = 'Sum';
        sum.type = 'Boolean';
        sum.default = 1;
    option.sum = sum;

        chwr.key = 'HalfwaveCenter';
        chwr.type = 'Boolean';
        chwr.default = 0;
        chwr.when = 'After';
    option.chwr = chwr;
    
        mu.key = 'Mu';
        mu.type = 'Integer';
        mu.default = 0;
        mu.keydefault = 100;
    option.mu = mu;
    
        oplog.key = 'Log';
        oplog.type = 'Boolean';
        oplog.default = 0;
        oplog.when = 'After';
    option.log = oplog;

        minlog.key = 'MinLog';
        minlog.type = 'Integer';
        minlog.default = 0;
        minlog.when = 'After';
    option.minlog = minlog;

        oppow.key = 'Power';
        oppow.type = 'Boolean';
        oppow.default = 0;
        oppow.when = 'After';
    option.power = oppow;
    
        diffenv.key = 'DiffEnvelope'; % obsolete, replaced by 'Diff'
        diffenv.type = 'Boolean';
        diffenv.default = 0;
    option.diffenv = diffenv;

        diff.key = 'Diff';
        diff.type = 'Integer';
        diff.default = 0;
        diff.keydefault = 1;
        diff.when = 'After';
    option.diff = diff;
    
        diffhwr.key = 'HalfwaveDiff';
        diffhwr.type = 'Integer';
        diffhwr.default = 0;
        diffhwr.keydefault = 1;
        diffhwr.when = 'After';
    option.diffhwr = diffhwr;

        lambda.key = 'Lambda';
        lambda.type = 'Integer';
        lambda.default = 1;
        lambda.when = 'After';
    option.lambda = lambda;

        c.key = 'Center';
        c.type = 'Boolean';
        c.default = 0;
        c.when = 'After';
    option.c = c;
    
        aver.key = 'Smooth';
        aver.type = 'Integer';
        aver.default = 0;
        aver.keydefault = 30;
        aver.when = 'After';
    option.aver = aver;
    
        ds.key = {'Down','PostDecim'};
        ds.type = 'Integer';
        if isamir(x,'mirenvelope')
            ds.default = 1;
        else
            ds.default = NaN;
        end
        ds.when = 'Both';
        ds.chunkcombine = 'During';
    option.ds = ds;

        sampling.key = 'Sampling';
        sampling.type = 'Integer';
        sampling.default = 0;
        sampling.when = 'After';
    option.sampling = sampling;
    
        up.key = {'UpSample'};
        up.type = 'Integer';
        up.default = 0;
        up.keydefault = 2;
    option.up = up;
    
        normal.key = 'Normal';
        normal.type = 'String';
        normal.choice = {0,1,'AcrossSegments'};
        normal.default = 1;
        normal.when = 'After';
    option.normal = normal;

%% options related to 'SpectralFlux'
        flux.key = 'SpectralFlux';
        flux.type = 'Boolean';
        flux.default = 0;
    option.flux = flux;
    
        complex.key = 'Complex';
        complex.type = 'Boolean';
        complex.when = 'Both';
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
        median.when = 'After';
    option.median = median;

        hw.key = 'Halfwave';
        hw.type = 'Boolean';
        hw.default = 1;
        hw.when = 'After';
    option.hw = hw;
    
%% options related to 'Pitch':
        pitch.key = 'Pitch';
        pitch.type = 'Boolean';
        pitch.default = 0;
    option.pitch = pitch;

        min.key = 'Min';
        min.type = 'Integer';
        min.default = 30;
    option.min = min;

        max.key = 'Max';
        max.type = 'Integer';
        max.default = 1000;
    option.max = max;
    
        novelty.key = 'Novelty';
        novelty.type = 'Boolean';
        novelty.default = 0;
    option.novelty = novelty;

        kernelsize.key = 'KernelSize';
        kernelsize.type = 'Integer';
        kernelsize.default = 0;
    option.kernelsize = kernelsize;
    
%% options related to 'Emerge':
        sgate.key = {'SmoothGate','Emerge'};
        sgate.type = 'String';
        sgate.choice = {'Goto','Lartillot'};
        sgate.default = '';
        sgate.keydefault = 'Lartillot';
        sgate.when = 'Both';
    option.sgate = sgate;
    
        minres.key = 'MinRes';
        minres.type = 'Integer';
        minres.default = 10;
    option.minres = minres;

%%
        nomodif.key = 'NoModif';
        nomodif.type = 'Boolean';
        nomodif.default = 0;
    option.nomodif = nomodif;

%% options related to event detection
        detect.key = 'Detect';
        detect.type = 'String';
        detect.choice = {'Peaks','Valleys',0,'no','off'};
        detect.default = 'Peaks';
        detect.keydefault = 'Peaks';
        detect.when = 'After';
    option.detect = detect;
    
        cthr.key = 'Contrast';
        cthr.type = 'Integer';
        cthr.default = NaN;
        cthr.when = 'After';
    option.cthr = cthr;

        thr.key = 'Threshold';
        thr.type = 'Integer';
        thr.default = 0;
        thr.when = 'After';
    option.thr = thr;
    
        single.key = 'Single';
        single.type = 'Boolean';
        single.default = 0;
        single.when = 'After';
    option.single = single;

        attack.key = {'Attack','Attacks'};
        attack.type = 'String';
        attack.choice = {'Derivate','Effort'};
        attack.default = 0;
        attack.keydefault = 'Derivate';
        attack.when = 'Both';
    option.attack = attack;
    
        alpha.key = 'Alpha';
        alpha.type = 'Integer';
        alpha.default = 3.75;
        alpha.when = 'After';
    option.alpha = alpha;  
    
        new.key = 'New';
        new.default = 0;
        new.when = 'After';
    option.new = new;
        
        decay.key = {'Decay','Decays','Release','Releases'};
        decay.type = 'Boolean';
        decay.default = 0;
        decay.when = 'Both';
    option.decay = decay;
    
%% preselection
        presel.choice = {'Scheirer','Klapuri99'};
        presel.type = 'String';
        presel.default = 0;
    option.presel = presel;

            
%% 'Frame' option
        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.when = 'After';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [3 .1];
    option.frame = frame;

specif.option = option;

specif.eachchunk = 'Normal';
specif.combinechunk = 'Concat';
specif.extensive = 1;

specif.title = 'Onset curve'; %used for miroptions

varargout = mirfunction(@mironsets,x,varargin,nargout,specif,@init,@main);


%% INIT

function [y type] = init(x,option)
if iscell(x)
    x = x{1};
end
if option.nomodif
    y = x;
    return
end
if ischar(option.presel)
    if strcmpi(option.presel,'Scheirer')
        option.filtertype = 'Scheirer';
        option.filter = 'HalfHann';
        option.envmeth = 'Filter';
    elseif strcmpi(option.presel,'Klapuri99')
        option.filtertype = 'Klapuri';
        option.filter = 'HalfHann';
        option.envmeth = 'Filter';
        option.decim = 180;
        option.mu = 100;
    end
end
if option.diffenv
    option.env = 1;
end
if isnan(option.env)
    if option.flux || option.pitch || option.novelty || ...
            (ischar(option.sgate) && ~isempty(option.sgate))
        option.env = 0;
    else
        option.env = 1;
    end
end
if ~option.kernelsize
    if option.pitch
        option.kernelsize = 32;
    elseif option.novelty
        option.kernelsize = 64;
    end
end
if isamir(x,'miraudio')
    y = [];
    if option.env
        if strcmpi(option.envmeth,'Filter')
            if isnan(option.filter)
                if ischar(option.attack) || option.decay
                    option.filter = 'Butter';
                else
                    option.filter = 'IIR';
                end
            end
            if isnan(option.hilb)
                if ischar(option.attack) || option.decay
                    option.hilb = 1;
                else
                    option.hilb = 0;
                end
            end
            if isnan(option.fb)
                if ischar(option.attack) || option.decay
                    option.fb = 0;
                else
                    option.fb = 40;
                end
            end
            
            if option.fb>1
                fb = mirfilterbank(x,option.filtertype,'NbChannels',option.fb);
            else
                fb = x;
            end
            y = mirenvelope(fb,'Filter','FilterType',option.filter,...
                'Hilbert',option.hilb,...
                'Tau',option.tau,'CutOff',option.cutoff,...
                'UpSample',option.up,...
                'PreDecim',option.decim,'PostDecim',[0 option.ds],...
                'Mu',option.mu,...
                'PreSilence',option.presilence,...
                'PostSilence',option.postsilence);
        else
            if isnan(option.specframe)
                if ischar(option.attack) || option.decay
                    option.specframe = [.03 .02];
                else
                    option.specframe = [.1 .1];
                end
            end
            y = mirenvelope(x,'Spectro',...
                'Frame',option.specframe(1),option.specframe(2),...
                'PowerSpectrum',option.powerspectrum,...
                'TimeSmooth',option.timesmooth,...
                'Terhardt',option.terhardt,...
                'PreSilence',option.presilence,...
                'PostSilence',option.postsilence);
        end

    end
    if option.flux
        z = mirflux(x,'Inc',option.inc,'Complex',option.complex); %,'Dist','City'); %%%%%%%%%%%%%%%%%???
        if isempty(y)
            y = z;
        else
            y = y+z;
        end
    end
    if option.pitch
        [unused ac] = mirpitch(x,'Frame','Min',option.min,'Max',option.max);
        z = mirnovelty(ac,'KernelSize',option.kernelsize);
        if isempty(y)
            y = z;
        else
            y = y+z;
        end
    elseif option.novelty
        s = mirspectrum(x,'max',1000,'Frame',.05,.2,'MinRes',3,'dB');
        %c = mircepstrum(x,'Frame',.05,.2);
        %[p ac] = mirpitch(x,'Frame');
        z = mirnovelty(s,'KernelSize',option.kernelsize,... 'Flux',...  
                      ...'Distance','Euclidean',...
                      'Similarity','oneminus');
        if isempty(y)
            y = z;
        else
            y = y+z;
        end
    elseif ischar(option.sgate) && ~isempty(option.sgate)
        if strcmpi(option.sgate,'Goto')
            x = miraudio(x,'Sampling',22050);
            y = mirspectrum(x,'Frame',.04644,.25);
        else
            y = mirspectrum(x,'Frame',.05,.2,....
                            'MinRes',option.minres,'dB','max',5000);
            if option.minres < 1 && isa(y,'mirdesign')
                y = set(y,'ChunkSizeFactor',get(x,'ChunkSizeFactor')*5); %20/option.minres);
            end
        end
        y = mirflux(y,'Inc','BackSmooth',option.sgate,'Dist','Gate');
    %% other ideas
        %y = mircepstrum(x,'min',50,'Hz','max',600,'Hz','Frame',.05,.2);
        %y = mirnovelty(y);%,'Width',1000);
        %y = mirsimatrix(y,'Width',1000); %'Distance','NewGate'
        %y = mirpeaks(y,'Contrast',.1,'Threshold',.3);
        %y = mirautocor(x,'Freq','max',5000,'Hz','Frame',.05,.2);
    end
elseif (option.pitch && not(isamir(x,'mirscalar'))) ...
        || isamir(x,'mirsimatrix')
    y = mirnovelty(x,'KernelSize',option.kernelsize);
elseif isamir(x,'mirscalar') || isamir(x,'mirenvelope') || ...
        (isamir(x,'mirspectrum') && ischar(option.sgate) && ~isempty(option.sgate))
    y = x;
else
    y = mirflux(x,'Inc',option.inc,'Complex',option.complex); %Not used...
end
if ischar(option.attack) || option.decay
    z = mironsets(x,option.envmeth,...
        'PreSilence',option.presilence,'PostSilence',option.postsilence);
    y = {y,z};
end
type = 'mirenvelope';


%% MAIN

function o = main(o,option,postoption)
if iscell(o)
    postoption.new = o{2};
    o = o{1};
end
if not(isempty(option)) && ischar(option.presel)
    if strcmpi(option.presel,'Scheirer')
        postoption.sampling = 200;
        postoption.diffhwr = 1;
        option.sum = 0;
        postoption.detect = 0;
    elseif strcmpi(option.presel,'Klapuri99')
        postoption.diffhwr = 1;
        option.sum = 0;
        postoption.ds = 0;
        o2 = o;
    end
end
if not(isempty(option)) && option.diffenv
    postoption.diff = 1;
end
if isa(o,'mirenvelope')
    if isfield(postoption,'sampling') && postoption.sampling
        o = mirenvelope(o,'Sampling',postoption.sampling);
    elseif isfield(postoption,'ds') 
        if isnan(postoption.ds)
            if option.decim || strcmpi(option.envmeth,'Spectro')
                postoption.ds = 0;
            else
                postoption.ds = 16;
            end
        end
        if postoption.ds
            o = mirenvelope(o,'Down',postoption.ds);
        end
    end
end
if isfield(postoption,'cthr')
    if isa(o,'mirenvelope')
        if postoption.power
            o = mirenvelope(o,'Power');
        end
        if postoption.diff
            o = mirenvelope(o,'Diff',postoption.diff,...
                              'Lambda',postoption.lambda,...
                              'Complex',postoption.complex);
        end
        if postoption.diffhwr
            o = mirenvelope(o,'HalfwaveDiff',postoption.diffhwr,...
                              'Lambda',postoption.lambda,...
                              'Complex',postoption.complex);
        end
        if postoption.aver
            o = mirenvelope(o,'Smooth',postoption.aver);
        end    
        if postoption.chwr
            o = mirenvelope(o,'HalfwaveCenter');
        end
    elseif isa(o,'mirscalar') && strcmp(get(o,'Title'),'Spectral flux') && ...
            ischar(postoption.sgate) && ~isempty(postoption.sgate)
        if postoption.median
            o = mirflux(o,'Median',postoption.median(1),postoption.median(2),...
                          'Halfwave',postoption.hw);
        else
            o = mirflux(o,'Halfwave',postoption.hw);
        end
    elseif isa(o,'mirscalar') && strcmp(get(o,'Title'),'Novelty')
        if postoption.diff
            o = mirenvelope(o,'Diff',postoption.diff,...
                              'Lambda',postoption.lambda,...
                              'Complex',postoption.complex);
        end
    end
end
if isa(o,'mirspectrum')
    [tmp o] = gettmp(o);
    d = get(o,'Data');
    [do tmp] = mircompute(@newonset,d,tmp);
    o = mirscalar(o,'Data',do,'Title','Onset curve');
    o = settmp(o,tmp);
%elseif isa(o,'mircepstrum')
%    pp = get(o,'PeakPosUnit');
%    pv = get(o,'PeakVal');
%    do = mircompute(@cepstronset,pp,pv);
%    o = mirscalar(o,'Data',do,'Title','Onset curve');
end

if isfield(option,'sum') && option.sum
    o = mirsum(o,'Adjacent',option.sum);
end
if isa(o,'mirenvelope') && isfield(postoption,'normal') && ...
        ~isequal(postoption.normal,0) && ~get(o,'Log')
    o = mirenvelope(o,'Normal',postoption.normal);
end
if isa(o,'mirenvelope') && isfield(postoption,'log') && postoption.log
    o = mirenvelope(o,'Log');
end
if isfield(option,'presel') && ...
        ischar(option.presel) && strcmpi(option.presel,'Klapuri99')
    % o, already computed, corresponds to mirenvelope(o,'Mu','HalfwaveDiff');
    % o is the relative distance function W in (Klapuri, 99);
    o2 = mirenvelope(o2,'HalfwaveDiff');
    % o2 is the absolute distance function D in (Klapuri, 99);
    p = mirpeaks(o,'Contrast',.2,'Chrono');
    p2 = mirpeaks(o2,'ScanForward',p,'Chrono');
    o = combinepeaks(p,p2,.05);
    clear o2 p p2
    filtfreq = 44*[2.^ ([ 0:2, ( 9+(0:17) )/3 ]) ];% Center frequencies of bands
    o = mirsum(o,'Weights',(filtfreq(1:end-1)+filtfreq(2:end))/2);
    o = mirenvelope(o,'Smooth',12);
end
if isfield(postoption,'detect')
    if postoption.c || (ischar(postoption.sgate) && ~isempty(postoption.sgate))
        o = mirenvelope(o,'Center');
    end
    if isa(o,'mirenvelope') && postoption.minlog
        o = mirenvelope(o,'MinLog',postoption.minlog);
    end
end
o = mirframenow(o,postoption);
if isfield(postoption,'detect') && ischar(postoption.detect)
    if isnan(postoption.cthr) || not(postoption.cthr)
        if ischar(postoption.attack) || postoption.decay
            postoption.cthr = .05;
        elseif ischar(postoption.detect) || postoption.detect
            postoption.cthr = .01;
        end
    elseif postoption.cthr
        if not(ischar(postoption.detect) || postoption.detect)
            postoption.detect = 'Peaks';
        end
    end
    if postoption.single
        total = 1;
        noend = 0;
    else
        total = Inf;
        noend = 1;
    end
    if strcmpi(postoption.detect,'Peaks')
        o = mirpeaks(o,'Total',total,'SelectFirst',0,...
            'Threshold',postoption.thr,'Contrast',postoption.cthr,...
            'Order','Abscissa','NoBegin','NoEnd',noend);
    elseif strcmpi(postoption.detect,'Valleys')
        o = mirpeaks(o,'Total',total,'SelectFirst',0,...
            'Threshold',postoption.thr,'Contrast',postoption.cthr,...
            'Valleys','Order','Abscissa','NoBegin','NoEnd',noend);
    end
    nop = cell(size(get(o,'Data')));
    o = set(o,'OnsetPos',nop,'AttackPos',nop,'DecayPos',nop);
end
if (isfield(postoption,'attack')) && (ischar(postoption.attack) || postoption.decay)
    pp = get(o,'PeakPos');
    d = get(o,'Data');
    t = get(o,'Time');
    if ischar(postoption.attack)
        x = postoption.new;
        ppu = get(o,'PeakPosUnit');
        if isnumeric(x)
            st = {{{}}};
            ap = {{{}}};
        else
            v = mirpeaks(x,'Total',Inf,'SelectFirst',0,...
                'Contrast',.1,...postoption.cthr,...
                'Threshold',.5,...
                'Valleys','Order','Abscissa','NoEnd');
            stu = get(v,'PeakPosUnit');
            [st,ap] = mircompute(@startattack,d,t,pp,ppu,stu,postoption);
        end
    else
        st = {{{}}};
        ap = {{{}}};
    end
    if postoption.decay
        x = postoption.new;
        ppu = get(o,'PeakPosUnit');
        if isnumeric(x)
            rl = {{{}}};
            en = {{{}}};
        else
            v = mirpeaks(x,'Total',Inf,'SelectFirst',0,...
                'Contrast',.1,...postoption.cthr,
                'Threshold',.5,...
                'Valleys','Order','Abscissa','NoBegin');
            rlu = get(v,'PeakPosUnit');
            [rl,en] = mircompute(@enddecay,d,t,pp,ppu,rlu);
        end
    else
        rl = {{{}}};
        en = {{{}}};
    end
    o = set(o,'OnsetPos',st,'AttackPos',ap,'DecayPos',rl,'OffsetPos',en,'PeakPos',pp);
end
title = get(o,'Title');
if not(length(title)>11 && strcmp(title(1:11),'Onset curve'))
    o = set(o,'Title',['Onset curve (',title,')']);
end


%%
function [do tmp] = newonset(d,tmp)
d = d - max(max(max(d)));
do = zeros(1,size(d,2));
if isempty(tmp)
    activ = [];
    inactiv = [];
    old = [];
else
    activ = tmp.activ;
    inactiv = tmp.inactiv;
    old = tmp.old;
end
for i = 1:size(d,2)
    dd = diff(d(:,i));
    new = find(dd(1:end-1) > 0 & dd(2:end) < 0) + 1;
    oldnew = [new d(new,i)];
    if ~isempty(old)
        maj = find(d(new,i) > -20);
        while ~isempty(maj)
            [min_o best_o] = min(abs(old(:,1) - new(maj(1))));
            min_a = Inf;
            for k = 1:length(activ)
                da = abs(activ(k).idx(end) - new(maj(1)));
                if da < min_a
                    min_a = da;
                    best_a = k;
                end
            end
            if min_a == min_o
                activ(best_a).idx(end+1) = new(maj(1));
                activ(best_a).mag(end+1) = d(new(maj(1)),i);
                activ(best_a).tim(end+1) = i;
                if length(activ(best_a).idx) < 10
                    do(i) = do(i) + d(new(maj(1)),i) + 20;
                end
            elseif old(best_o,2) < -20
                found = 0;
                for k = 1:length(inactiv)
                    if inactiv(k).idx(end) == old(best_o,1)
                        activ(end+1) = inactiv(k);
                        inactiv(k) = [];
                        activ(end).idx(end+1) = new(maj(1));
                        activ(end).mag(end+1) = d(new(maj(1)),i);
                        activ(end).tim(end+1) = i;
                        if length(activ(end).idx) < 10
                            do(i) = do(i) + d(new(maj(1)),i) + 20;
                        end
                        found = 1;
                        break
                    end
                end
                if ~found
                    activ(end+1).idx = new(maj(1));
                    activ(end).mag = d(new(maj(1)),i);
                    activ(end).tim = i;
                    do(i) = do(i) + d(new(maj(1)),i) + 20;
                end
            end
            new(maj(1)) = [];
            maj(1) = [];
            maj = maj - 1;
        end
                
        j = 1;
        while j <= length(inactiv)
            if isempty(new)
                inactiv(j:end) = [];
                break
            end
            
            [unused best_i] = min(abs(new - inactiv(j).idx(end)));
            if unused < 100
                inactiv(j).idx(end+1) = new(best_i);
                inactiv(j).mag(end+1) = d(new(best_i));
                inactiv(j).tim(end+1) = i;
                j = j+1;
            else
                inactiv(j) = [];
            end
            new(best_i) = [];
        end
        
        j = 1;
        while j <= length(activ)
            if activ(j).tim(end) < i
                if isempty(new)
                    activ(j) = [];
                    continue
                end
                if isempty(inactiv)
                    inactiv = activ(j);
                else
                    inactiv(end+1) = activ(j);
                end
                activ(j) = [];
                [unused best_i] = min(abs(new - inactiv(end).idx(end)));
                inactiv(end).idx(end+1) = new(best_i);
                inactiv(end).mag(end+1) = d(new(best_i));
                inactiv(end).tim(end+1) = i;
                new(best_i) = [];
            else
                j = j+1;
            end
        end
    end
    old = oldnew;
end
tmp.activ = activ;
tmp.inactiv = inactiv;
tmp.old = old;


%%
function [st pp] = startattack(d,t,pp,ppu,stu,option)
pp = sort(pp{1});
ppu = sort(ppu{1});
if isempty(pp)
    st = {{} {}};
    return
end

stu = stu{1};
if ~isempty(stu) && stu(1)>ppu(1)
    stu = [t(1) stu];
end
st = zeros(1,length(stu));

i = 0;
while i < length(stu)
    if length(ppu) == i
        break
    end
    i = i+1;
    
    % Removing additional peaks before current onset time stu(i)
    j = find(ppu(i:end) > stu(i),1);
    if j > 1
        ppu(i:i+j-2) = [];
        pp(i:i+j-2) = [];
    end
    
    % Taking the latest possible onset time before current peak ppu(i)
    j = find(stu(i:end) > ppu(i),1);
    if j > 2
        st(i:i+j-3) = [];
        stu(i:i+j-3) = [];
    end
    
    st(i) = find(t >= stu(i),1);
    
    while true
        dd = diff(d(st(i):pp(i)));
        f0 = find(dd > 0,1);
        if isempty(f0)
            pp(i) = [];
            ppu(i) = [];
        else
            break
        end
    end
    st(i) = st(i) + f0 - 1;
    
    if strcmpi(option.attack,'Derivate')
        dd = diff(d(st(i):pp(i)));
        f0 = find(dd < 0 & d(st(i):pp(i)-1) < d(st(i)));
        if ~isempty(f0)
            st(i) = st(i) + f0(end) - 1;
        end

        ppi = find(t >= ppu(i),1);
        dd = diff(d(st(i):pp(i)));
        f0 = find(dd <= 0 & ...
                  d(st(i):pp(i)-1) - d(st(i)) > (d(ppi) - d(st(i))) / 5 ...
                  ,1);
        if ~isempty(f0)
            pp(i) = st(i) + f0 - 1;
        end

        dd = diff(d(st(i):pp(i)));
        [mad, mdd] = max(dd);

        f2 = find(dd(end:-1:mdd+1)>mad/5,1);
        if isempty(f2)
            f2 = 1;
        end
        pp(i) = st(i) + length(dd) - f2;

        f1 = find(dd(1:mdd-1)>mad/10,1);
        if isempty(f1)
            f1 = 1;
        end
        st(i) = st(i) + f1;
    elseif strcmpi(option.attack,'Effort')
        % from Timbre Toolbox
        f_Env_v = d(st(i):pp(i));
        f_EnvMax = max(f_Env_v);
        f_Env_v = f_Env_v /f_EnvMax; % normalize by maximum value
        
        % === calcul de la pos pour chaque seuil
        percent_step	= 0.1;
        percent_value_v = percent_step:percent_step:1;
        percent_posn_v	= zeros(size(percent_value_v));
        for p=1:length(percent_value_v)
            percent_posn_v(p) = find(f_Env_v >= percent_value_v(p),1);
        end

        % === NOTATION
        % satt: start attack
        % eatt: end attack

        % === PARAMETRES
        param.m1	= 3; % === BORNES pour calcul mean
        param.m2	= 6;
        
        param.s1att	= 1; % === BORNES pour correction satt (start attack)
        param.s2att	= 3;
        
        param.e1att	= round(0.5/percent_step); % === BORNES pour correction eatt (end attack)
		param.e2att	= round(0.9/percent_step);
        
        % === dpercent_posn_v = effort
        dpercent_posn_v	= diff(percent_posn_v);
        % === M = effort moyen
        M				= mean(dpercent_posn_v(param.m1:param.m2));
        
        % === 1) START ATTACK
        % === on DEMARRE juste APRES que l'effort à fournir (écart temporal entre percent) soit trop important
        pos2_v			= find(dpercent_posn_v(param.s1att:param.s2att) > option.alpha * M);
        if ~isempty(pos2_v)
            result		= pos2_v(end)+param.s1att-1+1;
        else
            result		= param.s1att;
        end
        satt_posn		= percent_posn_v(result);
        
        % === raffinement: on cherche le minimum local
        delta	= round(0.25*(percent_posn_v(result+1)-percent_posn_v(result)));
        n		= percent_posn_v(result);
        if n-delta >= 1
            [min_value, min_pos]= min(f_Env_v(n-delta:n+delta));
            satt_posn			= min_pos + n-delta-1;
        end
        
        % === 2) END ATTACK
		% === on ARRETE juste AVANT que l'effort à fournir (écart temporal entre percent) soit trop important
		pos2_v		= find(dpercent_posn_v(param.e1att:param.e2att) > option.alpha * M);
		if ~isempty(pos2_v)
            result		= pos2_v(1)+param.e1att-1;
        else
            result		= param.e2att+1;
		end
		eatt_posn	= percent_posn_v(result);
        
		% === raffinement: on cherche le maximum local
		delta	= round(0.25*(percent_posn_v(result)-percent_posn_v(result-1)));
		n		= percent_posn_v(result);
		if n+delta <= length(f_Env_v)
			[max_value, max_pos]	= max(f_Env_v(n-delta:n+delta));
			eatt_posn				= max_pos + n-delta-1;
        end
        pp(i) = st(i) + eatt_posn;
        st(i) = st(i) + satt_posn;
    end
end
pp(length(st)+1:end) = [];
st = {{st} {pp}};


%%
function [pp en] = enddecay(d,t,pp,ppu,rlu)
pp = sort(pp{1});
ppu = sort(ppu{1});
if isempty(pp)
    en = {{} {}};
    return
end

rlu = rlu{1};
if ~isempty(rlu) && rlu(end)<ppu(end)
    rlu = [rlu t(end)];
end
en = zeros(1,length(rlu));

i = 0;
while i < length(rlu)
    if length(ppu) == i
        break
    end
    i = i+1;
    
     % Removing additional offset times before current peak ppu(i)
    j = find(rlu(i:end) > ppu(i),1);
    if j > 1
        rlu(i:i+j-2) = [];
        en(i:i+j-2) = [];
    end
    
    % Taking the latest possible peak before current offset time rlu(i)
    j = find(ppu(i:end) > rlu(i),1);
    if j > 2
        pp(i:i+j-3) = [];
        ppu(i:i+j-3) = [];
    end
        
    en(i) = find(t <= rlu(i),1,'last');
    
    stop = 0;
    while true
        if pp(i) >= en(i)
            stop = 1;
            break
        end
        dd = diff(d(en(i):-1:pp(i)));
        f0 = find(dd > 0,1);
        if isempty(f0)
            pp(i) = [];
            ppu(i) = [];
        else
            break
        end
    end
    if stop
        continue
    end
    en(i) = en(i) - f0 + 1;
    
    dd = diff(d(en(i):-1:pp(i)));
    f0 = find(dd < 0 & d(en(i):-1:pp(i)+1) < d(en(i)));
    if ~isempty(f0)
        en(i) = en(i) - f0(end) + 1;
    end
    
    ppi = find(t >= ppu(i),1);
    dd = diff(d(en(i):-1:pp(i)));
    f0 = find(dd <= 0 & ...
              d(en(i):-1:pp(i)+1) - d(en(i)) > (d(ppi) - d(en(i))) / 5 ...
              ,1);
    if ~isempty(f0)
        pp(i) = en(i) - f0 + 1;
    end
    
    dd = diff(d(en(i):-1:pp(i)));
    [mad, mdd] = max(dd);
    
    f2 = find(dd(end:-1:mdd+1)>mad/5,1);
    if isempty(f2)
        f2 = 1;
    end
    pp(i) = en(i) - length(dd) + f2;
    
    f1 = find(dd(1:mdd-1)>mad/10,1);
    if isempty(f1)
        f1 = 1;
    end
    en(i) = en(i) - f1;
end
pp(length(en)+1:end) = [];

pp = {{pp} {en}};