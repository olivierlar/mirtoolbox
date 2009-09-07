function varargout = mironsets(x,varargin)
%   o = mironsets(x) shows a temporal curve where peaks relate to the 
%       position of note onset times, and estimates those note onset 
%       positions.
%   Optional arguments:
%       mironsets(...,f) selects the strategy for the computation of the
%           onset detection function.
%           f = 'Envelope': Envelope of the audio signal. (Default choice).
%           With two methods for envelope extraction:
%               mironsets(...,'Spectro') (Default):
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
%                   then end of the computation. The resulting onset curve
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
%           f = 'Pitch ':computes a frame-decomposed autocorrelation function ,
%                of same default characteristics than those returned
%                 by mirpitch (with however a range of frequencies not
%                 exceeding 1000 Hz) and subsequently computes the novelty 
%                 curve of the resulting similatrix matrix, with a 'KernelSize' 
%                 of 32 samples.
%       mironsets(...,'Detect',d) toggles on or off the onset detection, 
%           which is based on the onset detection function.
%           (By default toggled on.)
%           Option associated to the mirpeaks function can be specified as
%               well:
%               'Contrast' with default value c = .01
%       mironsets(...,'Attack') (or 'Attacks') detects attack phases.
%       mironsets(...,'Release') (or 'Releases') detects release phases.
%           mironsets(...,'Gauss',o) estimate the attack and/or release
%               points using a gaussian envelope smoothing of order o of the
%               onset curve.
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
            filter.choice = {'IIR','HalfHann'};
            filter.default = 'IIR';
        option.filter = filter;

            tau.key = 'Tau';
            tau.type = 'Integer';
            tau.default = .02;
        option.tau = tau;

            fb.key = {'Filterbank','NbChannels'};
            fb.type = 'Integer';
            fb.default = 40;
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
    
%%      options related to 'Spectro':

            band.type = 'String';
            band.choice = {'Freq','Mel','Bark','Cents'};
            band.default = 'Freq';
        option.band = band;

                    
        sum.key = 'Sum';
        sum.type = 'Boolean';
        sum.default = 1;
    option.sum = sum;

        chwr.key = 'HalfwaveCenter';
        chwr.type = 'Boolean';
        chwr.default = 0;
        chwr.when = 'After';
    option.chwr = chwr;
    
        oplog.key = 'Log';
        oplog.type = 'Boolean';
        oplog.default = 0;
        oplog.when = 'After';
    option.log = oplog;
    
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
        ds.when = 'After';
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
    
        attack.key = {'Attack','Attacks'};
        attack.type = 'Boolean';
        attack.default = 0;
        attack.when = 'After';
    option.attack = attack;
        
        release.key = {'Release','Releases'};
        release.type = 'Boolean';
        release.default = 0;
        release.when = 'After';
    option.release = release;

        gauss.key = 'Gauss';
        gauss.type = 'Integer';
        gauss.default = 0;
        gauss.when = 'After';
    option.gauss = gauss;
    
%% preselection
        presel.choice = {'Scheirer','Klapuri99'};
        presel.type = 'String';
        presel.default = 0;
    option.presel = presel;

            
%% 'Frame' option
        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.when = 'Both';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [3 .1];
    option.frame = frame;

specif.option = option;

specif.eachchunk = 'Normal';
specif.combinechunk = 'Concat';

specif.title = 'Onset curve'; %used for miroptions
varargout = mirfunction(@mironsets,x,varargin,nargout,specif,@init,@main);


%% INIT

function [y type] = init(x,option)
if ischar(option.presel)
    if strcmpi(option.presel,'Scheirer')
        option.filtertype = 'Scheirer';
        option.filter = 'HalfHann';
    elseif strcmpi(option.presel,'Klapuri99')
        option.filtertype = 'Klapuri';
        option.filter = 'HalfHann';
        option.envmeth = 'Filter';
        option.decim = 180;
    end
end
if option.diffenv
    option.env = 1;
end
if isnan(option.env)
    if option.flux || option.pitch
        option.env = 0;
    else
        option.env = 1;
    end
end
if isamir(x,'miraudio')
    if option.env
        if strcmpi(option.envmeth,'Filter') && option.fb>1
            fb = mirfilterbank(x,option.filtertype,'NbChannels',option.fb);
        else
            fb = x;
        end
        y = mirenvelope(fb,option.envmeth,option.band,...
                          'FilterType',option.filter,...
                          'Tau',option.tau,'UpSample',option.up,...
                          'PreDecim',option.decim,'PostDecim',0);
        type = 'mirenvelope';
    elseif option.flux
        x = mirframenow(x,option);
        y = mirflux(x,'Inc',option.inc,'Complex',option.complex);
        type = 'mirscalar';
    elseif option.pitch
        [unused ac] = mirpitch(x,'Frame','Stable','Multi','Max',1000);
        y = mirnovelty(ac,'KernelSize',32);
        type = 'mirscalar';
    end
elseif (option.pitch && not(isamir(x,'mirscalar'))) ...
        || isamir(x,'mirsimatrix')
    y = mirnovelty(x,'KernelSize',32);
    type = 'mirscalar';
elseif isamir(x,'mirscalar') || isamir(x,'mirenvelope')
    y = x;
    type = mirtype(x);
else
    x = mirframenow(x,option);
    y = mirflux(x,'Inc',option.inc,'Complex',option.complex); %Not used...
    type = 'mirscalar';
end


%% MAIN

function o = main(o,option,postoption)
if not(isempty(option)) && ischar(option.presel)
    if strcmpi(option.presel,'Scheirer')
        postoption.sampling = 200;
        postoption.diffhwr = 1;
        option.sum = 0;
        postoption.detect = 0;
    elseif strcmpi(option.presel,'Klapuri99')
        postoption.log = 1;
        postoption.diffhwr = 1;
        option.sum = 0;
        postoption.ds = 0;
        o2 = o;
    end
end
if iscell(o)
    o = o{1};
end
if isfield(option,'sum') && option.sum && 0
    o = mirsum(o,'Adjacent',option.sum);
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
        if postoption.log
            o = mirenvelope(o,'Log');
        end
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
        if postoption.c
            o = mirenvelope(o,'Center');
        end
    elseif isa(o,'mirscalar') && strcmp(get(o,'Title'),'Spectral flux')
        if postoption.median
            o = mirflux(o,'Median',postoption.median(1),postoption.median(2),...
                          'Halfwave',postoption.hw);
        else
            o = mirflux(o,'Halfwave',postoption.hw);
        end
    end
end
if isfield(option,'sum') && option.sum
    o = mirsum(o,'Adjacent',option.sum);
end
if isfield(option,'presel') && ...
        ischar(option.presel) && strcmpi(option.presel,'Klapuri99')
    % o, already computed, corresponds to mirenvelope(o,'Log','HalfwaveDiff');
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
if not(isa(o,'mirscalar'))
    o = mirframenow(o,postoption);
end
if isfield(postoption,'detect') && ischar(postoption.detect)
    if isnan(postoption.cthr) || not(postoption.cthr)
        if ischar(postoption.detect) || postoption.detect
            postoption.cthr = .01;
        end
    elseif postoption.cthr
        if not(ischar(postoption.detect) || postoption.detect)
            postoption.detect = 'Peaks';
        end
    end
    if strcmpi(postoption.detect,'Peaks')
        o = mirpeaks(o,'Total',Inf,'SelectFirst',...
            'Contrast',postoption.cthr,'Order','Abscissa');
    elseif strcmpi(postoption.detect,'Valleys')
        o = mirpeaks(o,'Total',Inf,'SelectFirst',...
            'Contrast',postoption.cthr,'Valleys','Order','Abscissa');
    end
    nop = cell(size(get(o,'Data')));
    o = set(o,'AttackPos',nop,'ReleasePos',nop);
end
if isfield(postoption,'attack') && postoption.attack
    p = get(o,'PeakPos');
    pm = get(o,'PeakMode');
    d = get(o,'Data');
    do = mirenvelope(o,'Diff','Gauss',postoption.gauss);
    dd = get(do,'Data');
    [st p pm] = mircompute(@startattack,d,dd,p,pm);
    o = set(o,'AttackPos',st,'PeakPos',p,'PeakMode',pm);
else
    st = [];
end
if isfield(postoption,'release') && postoption.release
    p = get(o,'PeakPos');
    pm = get(o,'PeakMode');
    d = get(o,'Data');
    do = mirenvelope(o,'Diff','Gauss',postoption.gauss);
    dd = get(do,'Data');
    [rl p pm st] = mircompute(@endrelease,d,dd,p,pm,st);
    o = set(o,'ReleasePos',rl,'AttackPos',st,'PeakPos',p,'PeakMode',pm);
end
title = get(o,'Title');
if not(length(title)>11 && strcmp(title(1:11),'Onset curve'))
    o = set(o,'Title',['Onset curve (',title,')']);
end


function st = startattack(d,dd,z,pm)
z = sort(z{1});
pm = pm{1};
st = zeros(size(z));
i = 1;
while i<=length(z)
    p = find(dd(z(i)-1:-1:1)>0,1); % previous increasing portion
    if isempty(p)
        st(i) = 1;
    else
        n = find(dd(z(i)-p-1:-1:1)<0,1); % previous decreasing portion
        if isempty(n)
            st(i) = 1;
        else
            st(i) = z(i)-p-n;
        end
        if i>1 && st(i-1)==st(i)
            if d(z(i))>d(z(i-1))
                del = i-1;
            else
                del = i;
            end
            st(del) = [];
            z(del) = [];
            pm(del) = [];
            i = i-1;
        end
    end
    i = i+1;
end
st = {{st} {z} {pm}};


function rt = endrelease(d,dd,z,pm,st)
z = sort(z{1});
pm = pm{1};
if not(isempty(st))
    st = st{1};
end
rt = zeros(size(z));
i = 1;
while i<=length(z)
    p = find(dd(z(i)+1:end)<0,1); % next decreasing portion
    if isempty(p)
        rt(i) = length(d);
    else
        n = find(dd(z(i)+p+1:end)>0,1); % next increasing portion
        if isempty(n)
            rt(i) = length(d);
        else
            rt(i) = z(i)+p+n;
        end
        if i>1 && rt(i-1)==rt(i)
            if d(z(i))>d(z(i-1))
                del = i-1;
            else
                del = i;
            end
            rt(del) = [];
            z(del) = [];
            pm(del) = [];
            if not(isempty(st))
                st(del) = [];
            end
            i = i-1;
        end
    end
    i = i+1;
end
rt = {{rt} {z} {pm} {st}};