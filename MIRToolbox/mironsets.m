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
%       mironsets(...,'Release') (or 'Releases') detects release phases.
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
        
            specframe.key = 'SpectroFrame';
            specframe.type = 'Integer';
            specframe.number = 2;
            specframe.default = [.1 .1];
        option.specframe = specframe;
                    
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
        attack.type = 'Boolean';
        attack.default = 0;
        attack.when = 'After';
    option.attack = attack;
        
        release.key = {'Release','Releases'};
        release.type = 'String';
        release.choice = {'Olivier','Valeri',0,'no','off'};
        release.default = 0;
        release.keydefault = 'Olivier';
        release.when = 'After';
    option.release = release;
    
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
    if option.flux || option.pitch || option.novelty
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
        if strcmpi(option.envmeth,'Filter') && option.fb>1
            fb = mirfilterbank(x,option.filtertype,'NbChannels',option.fb);
        else
            fb = x;
        end
        y = mirenvelope(fb,option.envmeth,option.band,...
                          'Frame',option.specframe(1),option.specframe(2),...
                          'FilterType',option.filter,...
                          'Tau',option.tau,'UpSample',option.up,...
                          'PreDecim',option.decim,'PostDecim',0,...
                          'Mu',option.mu);
        type = 'mirenvelope';
    end
    if option.flux
        z = mirflux(x,'Inc',option.inc,'Complex',option.complex);
        if isempty(y)
            y = z;
        else
            y = y+z;
        end
        type = 'mirscalar';
    end
    if option.pitch
        [unused ac] = mirpitch(x,'Frame','Min',option.min,'Max',option.max);
        z = mirnovelty(ac,'KernelSize',option.kernelsize);
        if isempty(y)
            y = z;
        else
            y = y+z;
        end
        type = 'mirscalar';
    elseif option.novelty
        s = mirspectrum(x,'max',1000,'Frame',.05,.2,'MinRes',3);%,'log');
        %c = mircepstrum(x,'Frame',.05,.2);
        %[p ac] = mirpitch(x,'Frame');
        z = mirnovelty(s,'KernelSize',option.kernelsize,...
                      ...'Distance','Euclidean',...
                      'Similarity','oneminus');
        if isempty(y)
            y = z;
        else
            y = y+z;
        end
        type = 'mirscalar';
    end
elseif (option.pitch && not(isamir(x,'mirscalar'))) ...
        || isamir(x,'mirsimatrix')
    y = mirnovelty(x,'KernelSize',option.kernelsize);
    type = 'mirscalar';
elseif isamir(x,'mirscalar') || isamir(x,'mirenvelope')
    y = x; %mirframenow(x,option);
    type = mirtype(x);
else
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
        postoption.diffhwr = 1;
        option.sum = 0;
        postoption.ds = 0;
        o2 = o;
    end
end
if iscell(o)
    o = o{1};
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
    o = mirenvelope(o,'Normal');
end
o = mirframenow(o,postoption);
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
    o = set(o,'AttackPos',nop,'ReleasePos',nop);
end
if (isfield(postoption,'attack') && not(isequal(postoption.attack,0))) || ...
        (isfield(postoption,'release') && not(isequal(postoption.release,0)))
    pp = get(o,'PeakPos');
    pv = get(o,'PeakVal');
    pm = get(o,'PeakMode');
    ppp = get(o,'PeakPrecisePos');
    ppv = get(o,'PeakPreciseVal');
    d = get(o,'Data');
    if postoption.attack
        v = mirpeaks(o,'Total',Inf,'SelectFirst',0,...
            'Contrast',postoption.cthr,...
            'Valleys','Order','Abscissa','NoEnd');
        st = get(v,'PeakPos');
        [st pp] = mircompute(@startattack,d,pp,st);
        %[st pp pv pm ppp ppv] = mircompute(@startattack,d,pp,pv,pm,ppp,ppv);
    else
        st = [];
    end
    if ischar(postoption.release) && ~strcmpi(postoption.release,'No') ...
                                  && ~strcmpi(postoption.release,'Off')
        [rl pp pv pm ppp ppv st] = mircompute(@endrelease,d,pp,pv,pm,ppp,ppv,st,postoption.release);
        o = set(o,'ReleasePos',rl);
    end
    o = set(o,'AttackPos',st,'PeakPos',pp,'PeakVal',pv,'PeakMode',pm,...
              'PeakPrecisePos',ppp,'PeakPreciseVal',ppv);
end
title = get(o,'Title');
if not(length(title)>11 && strcmp(title(1:11),'Onset curve'))
    o = set(o,'Title',['Onset curve (',title,')']);
end


function st = startattack(d,pp,st) %pv,pm,ppp,ppv)
pp = sort(pp{1});
if isempty(pp)
    st = {{} {}};
    return
end

st = st{1};
if ~isempty(st) && st(1)>pp(1)
    dd = diff(d,1,1);       % d'
    p = find(dd((pp(1)-1)-1:-1:1)<=0,1);
    if isempty(p)
        st0 = 1;
    else
        st0 = ((pp(1)-1)-p)+1;
    end
    st = [st0 st];
end
if length(st)>length(pp)
    st = st(1:length(pp));
else
    pp = pp(1:length(st));
end
st = {{st} {pp}};
return

pv = pv{1};
pm = pm{1};
ppp = ppp{1};
ppv = ppv{1};
st = zeros(size(pp));
i = 1;
dd = diff(d,1,1);       % d'
ddd = diff(dd,1,1);     % d''
dddd = diff(ddd,1,1);   % d'''
while i<=length(pp)
    % Start attack is identified to previous valley in d.
    p = find(dd((pp(i)-1)-1:-1:1)<0,1);
    if isempty(p)
        st(i) = 1;
    else
        st(i) = ((pp(i)-1)-p)+1;
        if i>1 && st(i-1)==st(i)
            if d(pp(i))>d(pp(i-1))
                del = i-1;
            else
                del = i;
            end
            st(del) = [];
            pp(del) = [];
            pv(del) = [];
            pm(del) = [];
            ppp(del) = [];
            ppv(del) = [];
            i = i-1;
        end
    end
    % Start attack is identified to previous peak in d''.
    %p = find(dddd((pp(i)-1)-1:-1:1)<0,1); % previous decreasing d''
    %if isempty(p)
    %    st(i) = 1;
    %else
    %    n = find(dddd((pp(i)-1)-p-1:-1:1)>0,1); % previous increasing d''
    %    if isempty(n)
    %        st(i) = 1;
    %    else
    %        st(i) = ((pp(i)-1)-p-(n-1))+1;
    %    end
    %    if i>1 && st(i-1)==st(i)
    %        if d(pp(i))>d(pp(i-1))
    %            del = i-1;
    %        else
    %            del = i;
    %        end
    %        st(del) = [];
    %        pp(del) = [];
    %        pv(del) = [];
    %        pm(del) = [];
    %        ppp(del) = [];
    %        ppv(del) = [];
    %        i = i-1;
    %    end
    %end
    i = i+1;
end
st = {{st} {pp} {pv} {pm} {ppp} {ppv}};


function rt = endrelease(d,pp,pv,pm,ppp,ppv,st,meth)
if isempty(pp)
    rt = [];
    return
end
pp = sort(pp{1});
pv = pv{1};
pm = pm{1};
ppp = ppp{1};
ppv = ppv{1};
if not(isempty(st))
    st = st{1};
end
rt = zeros(size(pp));
i = 1;
dd = diff(d,1,1);       % d'
ddd = diff(dd,1,1);     % d''
dddd = diff(ddd,1,1);   % d'''
while i<=length(pp)
    if strcmpi(meth,'Olivier')
        % Release attack is identified to next (sufficiently positive) peak 
        % in d''.
        l = find(ddd((pp(i)-1):end)<min(ddd)/100,1); 
            % next d'' sufficiently negative
        if isempty(l)
            rt(i) = length(d);
        else
            p = find(ddd((pp(i)-1)+(l-1)+1:end)>max(ddd)/100,1); % next increasing d''
            if isempty(p)
                rt(i) = length(d);
            else
                n = find(dddd((pp(i)-1)+(l-1)+p+1:end)<0,1); % next decreasing d''
                if isempty(n)
                    rt(i) = length(d);
                else
                    rt(i) = ((pp(i)-1)+(l-1)+p+n)+1;
                end
            end
        end
    elseif strcmpi(meth,'Valeri')
        p = find(dd((pp(i)-1)+1:end)>min(dd)/100,1); % find point nearest to min(dd)/100 from current peak. 
        if isempty(p)
            rt(i) = length(d);
        elseif p<=3                                 %that means if p is less than 3 points away from the peak then it can not be considered as the end point of release.
                                                  %Assumption is that the whole DSR(decay sustain release) section can not be shorter than 30 ms (sampling rate is 100 Hz), also, no successive note can be nearer than 30ms.
            rt(i) = pp(i)+3;
        else
            rt(i) = (pp(i)-1)+(p-1);
        end
    end
    if i>1 && rt(i-1)==rt(i)
        if d(pp(i))>d(pp(i-1))
            del = i-1;
        else
            del = i;
        end
        rt(del) = [];
        pp(del) = [];
        pv(del) = [];
        pm(del) = [];
        ppp(del) = [];
        ppv(del) = [];
        if not(isempty(st))
            st(del) = [];
        end
        i = i-1;
    end
    i = i+1;
end
rt = {{rt} {pp} {pv} {pm} {ppp} {ppv} {st}};