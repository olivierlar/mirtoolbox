function varargout = mirmetre(orig,varargin)
%   m = mirmetre(x) estimates the metrical hierarchy.
%   Optional arguments:
%       mirtempo(...,'Frame',l,h) orders a frame decomposition of window
%           length l (in seconds) and hop factor h, expressed relatively to
%           the window length. For instance h = 1 indicates no overlap.
%           Default values: l = 3 seconds and h = .1
%       The options related to the onset detection phase can be specified 
%               (see help mironsets):
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
%
%   [m,p] = mirmetre(...) also displays the autocorrelation function
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
        frame.keydefault = [5 .05];
    option.frame = frame;
    
        fea.type = 'String';
        fea.choice = {'Envelope','DiffEnvelope','SpectralFlux',...
                      'Pitch','Novelty'};
        fea.default = 'Envelope'; %'Novelty';
    option.fea = fea;
    
    %% options related to 'Envelope':
    
            envmeth.key = 'Method';
            envmeth.type = 'String';
            envmeth.choice = {'Filter','Spectro'};
            envmeth.default = 'Spectro';
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
            diff.default = 1;
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
        
        
%% option related to mirautocor:                
    
        nw.key = 'NormalWindow';
        nw.default = 0;
    option.nw = nw;

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
        
        thr.key = 'Threshold';
        thr.type = 'Integer';
        thr.default = 0;
    option.thr = thr;
    
        cthr.key = 'Contrast';
        cthr.type = 'Integer';
        cthr.default = .05;
    option.cthr = cthr;

        mi.key = 'Min';
        mi.type = 'Integer';
        mi.default = 24;
    option.mi = mi;
        
        ma.key = 'Max';
        ma.type = 'Integer';
        ma.default = Inf;
    option.ma = ma;
    
        lart.key = 'Lartillot';
        lart.type = 'Integer';
        lart.default = .15;
    option.lart = lart;

        lart2.type = 'Integer';
        lart2.default = .45;
    option.lart2 = lart2;
        
specif.option = option;

varargout = mirfunction(@mirmetre,orig,varargin,nargout,specif,@init,@main);


function [y type] = init(x,option)
if iscell(x)
    x = x{1};
end
if isamir(x,'mirmetre')
    y = x;
    return
end
if not(isamir(x,'mirautocor')) && not(isamir(x,'mirspectrum'))
    if isframed(x) && strcmpi(option.fea,'Envelope') && not(isamir(x,'mirscalar'))
        warning('WARNING IN MIRMETRE: The input should not be already decomposed into frames.');
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
y = mirautocor(x,'Min',60/option.ma,'Max',60/option.mi * 2,...
      'NormalWindow',option.nw);
if ischar(option.sum)
    y = mirsum(y);
end
y = mirpeaks(y,'Total',Inf,...
               'Threshold',option.thr,'Contrast',option.cthr,...
               'NoBegin','NoEnd',...
               'Normalize','Local','Order','Amplitude');
type = {'mirmetre',mirtype(y)};
    

function o = main(p,option,postoption)
if iscell(p)
    p = p{1};
end
if isamir(p,'mirscalar')
    m = modif(m,postoption);
    o = {t};
    return
end
pt = get(p,'PeakPrecisePos');
meters = cell(1,length(pt));
d = get(p,'Data');
pp = get(p,'Pos');
ppp = get(p,'PeakPos');
for j = 1:length(pt)
    for k = 1:length(pt{j})
        ptk = pt{j}{k};
        for h = 1:size(ptk,3)
            mk = {};
            best = [];
            for l = 1:size(ptk,2)       % For each successive frame
                if 0 %l > 1
                    l
                    [mk{1}.lastbpm].*repmat([mk{1}.lvl],[2 1])
                end
                ptl = getbpm(p,ptk{1,l,h}); % Peaks
                cntr = zeros(1,length(mk));
                
                bpms = cell(1,length(mk));
                for i2 = 1:length(mk)
                    bpms{i2} = [mk{i2}.lastbpm];
                end
                ampli = d{j}{k}(:,l,h);
                pos = ppp{j}{k}{l};
                
                for i2 = 1:length(mk)   % For each metrical hierarchy
                    best(i2) = findbest(mk{i2},l-1);
                end
                
                for i = 1:length(ptl)       % For each peak
                    if ptl(i) < option.mi
                        continue
                    end
                      
                    if ~find(pos > pos(i)*1.95 & pos < pos(i)*2.05)
                        continue
                    end
                    
                    ptli = ptl(i);
                    delta1 = find(ampli(pos(i)+1:end) < 0,1);
                    if isempty(delta1)
                        delta1 = length(ampli) - pos(i);
                    end
                    delta2 = find(ampli(pos(i)-1:-1:1) < 0,1);
                    if isempty(delta2)
                        delta2 = pos(i) - 1;
                    end
                    ptli1 = getbpm(p,pp{j}{k}(pos(i)+delta1,l));
                    ptli2 = getbpm(p,pp{j}{k}(pos(i)-delta2,l));
                    
                    score = ampli(pos(i));
                    found = 0;              % Is peak in metrical hierarchies?
                    coord = [];             % Where is peak located
                    for i2 = 1:length(mk)   % For each metrical hierarchy
                        if ~isempty(bpms{i2})
                            locoord = [];
                            
                            bpm2 = repmat(bpms{i2}(:,best(i2))*mk{i2}(best(i2)).lvl,...
                                          [1 length(mk{i2})])...
                                   ./ repmat([mk{i2}.lvl],[2 1]);
                            for i3 = 1:size(bpm2,2)
                                if mk{i2}(i3).timidx < l-5
                                    bpm2(:,i3) = nan(2,1);
                                end
                            end
                            
                            dif = 60/ptli - 60./bpm2;
                            dif(:,( dif(1,:).*dif(2,:) < 0)) = 0;
                            dist = min(abs(dif));
                            [unused i3] = min(dist);
                            if ptli > bpm2(1,i3)
                                ratio = ptli/bpm2(1,i3);
                            else
                                ratio = bpm2(1,i3)/ptli;
                            end
                            if min(abs(60/ptli - 60/mk{i2}(i3).lastbpm(1))) < option.lart && ...
                                    mod(ratio,1) < option.lart2
                                % Continuing an existing metrical level.
                                
                                if mk{i2}(i3).timidx(end) == l
                                    % Already continued.
                                    if 0
                                        mk{i2}(end+1) = mk{i2}(i3);
                                        mk{i2}(end).bpms(end) = ptli;
                                        mk{i2}(end).score(end) = ...
                                            d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                        bpms{i2}(:,end+1) = [ptli; ptli]; %[ptli1; ptli2];
                                        locoord = length(mk{i2});
                                    end
                                    found = 1;
                                    
                                else
                                    if mean(bpms{i2}(:,i3)) > ptli
                                        div = mean(bpms{i2}(:,i3)) / ptl(i);
                                    else
                                        div = ptli / mean(bpms{i2}(:,i3));
                                    end
                                    % Continuing an existing metrical
                                    % level.
                                    found = 1;
                                    coord = [i2 i3];

                                    if isempty(locoord) || ...
                                            dist(i3) < dist(locoord)
                                        % Level not identified to
                                        % one already detected
                                        mk{i2}(i3).timidx(end+1) = l;
                                        mk{i2}(i3).bpms(end+1) = ptl(i);
                                        mk{i2}(i3).lastbpm = [ptli; ptli]; %[ptli1; ptli2];
                                        mk{i2}(i3).score(end+1) = ...
                                            d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                        %bpms{i2}(:,i3) = [ptli1; ptli2];
                                        if i == 1
                                            best(i2) = i3;
                                        end

                                        locoord = i3;
                                    end
                                end
                            end
                        end
                    end

                    if found
                        continue
                    end
                    
                    %%
                    % Candidate level not belonging to current
                    % metrical levels
                                        
                    incoherent = zeros(1,length(mk));
                    i2 = 1;
                    while i2 <= length(mk)
                        [unused ord] = sort(mean(bpms{i2}));
                        orbpms = bpms{i2}(:,ord);
                        orbest = find(ord == best(i2));

                        fo = find(orbpms(2,:) > ptli, 1);
                        if isempty(fo)
                            fo = size(orbpms,2)+1;
                        end

                        % Stored levels slower than candidate
                        slower = [];
                        i3 = fo-1;
                        err = Inf;
                        while i3 > 0
                            bpm3 = bpms{i2}(:,best(i2)) ...
                                        * mk{i2}(best(i2)).lvl ...
                                        / mk{i2}(ord(i3)).lvl;
                            div = [ptli2; ptli1] ./ bpm3; %orbpms(:,i3);
                            rdiv = round(ptli / mean(bpm3));
                            if rdiv > 1
                                if floor(div(1)) ~= floor(div(2))
                                    newerr = 0;
                                else
                                    newerr = min(min(mod(div,1)),...
                                             min(1-mod(div,1)));
                                end
                                if newerr < option.lart2
                                    % Candidate level can be
                                    % integrated in this metrical
                                    % hierarchy
                                    if newerr < err
                                        if isempty(slower)
                                            slower.lvl = mk{i2}(ord(i3)).lvl / rdiv;
                                            slower.bpm = mean(orbpms(:,i3));
                                            slower.score = mk{i2}(ord(i3)).score(end);
                                            slower.rdiv = rdiv;
                                            rptli1 = orbpms(1,i3) * (rdiv - .4);
                                            if ptli1 < rptli1
                                                ptli1 = rptli1;
                                            end
                                            rptli2 = orbpms(2,i3) * (rdiv + .4);
                                            if ptli2 > rptli2
                                                ptli2 = rptli2;
                                            end
                                            ptli = mean([ptli1,ptli2]);
                                            err = newerr;
                                            %break
                                        elseif mk{i2}(ord(i3)).lvl / rdiv ...
                                                ~= slower.lvl
                                            slower.lvl = mk{i2}(ord(i3)).lvl / rdiv;
                                            slower.rdiv = rdiv;
                                            %slower = [];
                                            %incoherent(i2) = 1;
                                            %break
                                        end
                                    end
                                end
                            elseif ~isempty(slower) && ...
                                    ~mod(mk{i2}(ord(i3)).lvl,slower.lvl)
                                slower = [];
                                incoherent(i2) = 1;
                                break
                            end
                            i3 = i3 - 1;
                        end

                        if incoherent(i2)
                            break
                        end

                        % Stored levels faster than candidate
                        faster = [];
                        i3 = fo;
                        err = Inf;
                        while i3 <= size(orbpms,2)
                            bpm3 = bpms{i2}(:,best(i2)) ...
                                        * mk{i2}(best(i2)).lvl ...
                                        / mk{i2}(ord(i3)).lvl;
                            div = bpm3 ./ [ptli2;ptli1];
                            rdiv = round(mean(bpm3) / ptli);
                            if rdiv > 1
                                if floor(div(1)) < floor(div(2))
                                    newerr = 0;
                                else
                                    newerr = min(min(mod(div,1)),...
                                             min(1-mod(div,1)));
                                end
                                if newerr < option.lart2
                                    % Candidate level can be
                                    % integrated in this metrical
                                    % hierarchy
                                    if newerr < err
                                        if isempty(faster)
                                            faster.lvl = mk{i2}(ord(i3)).lvl * rdiv;
                                            faster.bpm = mean(orbpms(:,i3));
                                            faster.score = mk{i2}(ord(i3)).score(end);
                                            faster.rdiv = rdiv;
                                            rptli1 = orbpms(1,i3) / (rdiv + .4);
                                            if ptli1 < rptli1
                                                ptli1 = rptli1;
                                            end
                                            rptli2 = orbpms(2,i3) / (rdiv - .4);
                                            if ptli2 > rptli2
                                                ptli2 = rptli2;
                                            end
                                            ptli = mean([ptli1,ptli2]);
                                            err = newerr;
                                            %break
                                        elseif mk{i2}(ord(i3)).lvl * rdiv ...
                                                ~= faster.lvl
                                            faster.lvl = mk{i2}(ord(i3)).lvl * rdiv;
                                            faster.rdiv = rdiv;
                                            %faster = [];
                                            %incoherent(i2) = 1;
                                            %break
                                        end
                                    end
                                end
                            elseif ~isempty(faster) && ...
                                    ~mod(faster.lvl,mk{i2}(ord(i3)).lvl)
                                faster = [];
                                incoherent(i2) = 1;
                                break
                            end
                            i3 = i3 + 1;
                        end

                        if incoherent(i2)
                            break
                        end

                        if isempty(slower) && isempty(faster)
                            i2 = i2 + 1;
                            continue
                        elseif isempty(slower)
                            lvl = faster.lvl;
                            rdiv = faster.rdiv;
                        elseif isempty(faster)
                            lvl = slower.lvl;
                            rdiv = slower.rdiv;
                        elseif slower.score < faster.score
                            lvl = faster.lvl;
                            rdiv = faster.rdiv;
                        else
                            lvl = slower.lvl;
                            rdiv = slower.rdiv;
                        end
                                                    
                        if ~found
                            found = 1;
                            l0 = find(lvl == [mk{i2}.lvl]);
                            if isempty(l0)
                                % New metrical level
                                mk{i2}(end+1).lvl = lvl;
                                mk{i2}(end).lastbpm = [ptli; ptli]; %[ptli1; ptli2];
                                mk{i2}(end).bpms = ptl(i);
                                mk{i2}(end).timidx = l;
                                mk{i2}(end).score = ampli(pos(i));
                                coord = [i2 length(mk{i2})];
                                bpms{i2}(:,end+1) = [ptli; ptli]; %[ptli1; ptli2];
                                if i == 1
                                    best(i2) = length(mk{i2});
                                end
                                if 0 %lvl<1
                                    for i4 = 1:length(mk{i2})
                                        mk{i2}(i4).lvl = ...
                                            mk{i2}(i4).lvl * rdiv;
                                    end
                                end
                            else
                                dist = abs(mean([mk{i2}(l0).lastbpm]) - ptl(i));
                                [unused md] = min(dist);
                                if mk{i2}(l0(md)).timidx(end) < l
                                    mk{i2}(l0(md)).lastbpm = [ptli; ptli]; %[ptli1; ptli2];
                                    mk{i2}(l0(md)).bpms(end+1) = ptl(i);
                                    mk{i2}(l0(md)).score(end+1) = ampli(pos(i));
                                    mk{i2}(l0(md)).timidx(end+1) = l;
                                elseif score > mk{i2}(l0(md)).score
                                    mk{i2}(l0(md)).lastbpm = [ptli; ptli]; %[ptli1; ptli2];
                                    mk{i2}(l0(md)).bpms(end) = ptl(i);
                                    mk{i2}(l0(md)).score(end) = ampli(pos(i));
                                end
                                coord = [i2 l0(md)];
                                if i == 1
                                    best(i2) = l0(md);
                                end
                            end
                            locoord = coord(2);
                            
                        else
                            % Candidate level also integrated in other
                            % metrical hierarchy. Both hierarchies are fused.
                            if lvl<1
                                for i4 = 1:length(mk{i2})
                                    mk{i2}(i4).lvl = mk{i2}(i4).lvl * rdiv;
                                end
                                lvl = lvl * rdiv;
                            end
                            
                            if mk{coord(1)}(coord(2)).lvl > lvl
                                % Other hierarchy is
                                % faster than current.
                                meter1 = mk{i2};
                                    % Slower hierarchy,
                                    % which is fused to
                                    % faster one
                                meter2 = mk{coord(1)};
                                    % Faster hierarchy,
                                    % onto which slower
                                    % one is fused
                                bpms2 = bpms{coord(1)};
                                lvl1 = lvl;
                                    % Level in slower
                                lvl2 = meter2(coord(2)).lvl;
                                    % Level in faster
                            else
                                % Other hierarchy is
                                % slower than current
                                meter1 = mk{coord(1)};
                                meter2 = mk{i2};
                                bpms2 = bpms{i2};
                                lvl1 = meter1(coord(2)).lvl;
                                lvl2 = lvl;
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
                                for i4 = 1:length(meter2)
                                    meter2(i4).lvl = ...
                                            meter2(i4).lvl * mult2;
                                end
                            end
                            for i4 = 1:length(meter1)
                                % Fusing one hierarchy
                                % into the other..
                                f = find(meter1(i4).lvl * mult1 ...
                                         == [meter2.lvl]);
                                if isempty(f)
                                    meter2(end+1) = meter1(i4);
                                    meter2(end).lvl = ...
                                        meter2(end).lvl * mult1;
                                    f = length(meter2);
                                    bpms2(:,end+1) = ...
                                                meter1(i4).lastbpm;
                                end
                                if coord(2) == i4
                                    coord(2) = f;
                                end
                            end
                            mk{coord(1)} = meter2;
                            bpms{coord(1)} = bpms2;
                            mk(i2) = [];
                            best(coord(1)) = findbest(meter2,l);
                            bpms(i2) = [];
                            best(i2) = [];
                            incoherent(i2) = [];
                            i2 = i2 - 1;
                        end
                        
                        i2 = i2 + 1;
                    end
                    
                    if ~found && isempty(find(incoherent))
                        % New metrical hierarchy
                        mk{end+1}.lvl = 1;
                        mk{end}.lastbpm = [ptli; ptli]; %[ptli1; ptli2];
                        mk{end}.bpms = ptl(i);
                        mk{end}.timidx = l;
                        mk{end}.score = ...
                            d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                        mk{end}.globpms = [];
                        cntr(end+1) = 0;
                        %found = 1;
                        bpms{end+1} = [ptli; ptli]; %[ptli1;ptli2];
                        best(end+1) = 1;
                        %coord = [length(bpms),1];
                    end
                end
                
                for i = 1:length(mk)
                    for i2 = 1:length(mk{i})
                        if mk{i}(i2).timidx(end) == l
                            mk{i}(i2).globpms(end+1) =  ...
                                mk{i}(best(i)).bpms(end) ...
                                    * mk{i}(best(i)).lvl ...
                                    / mk{i}(i2).lvl;
                        end
                    end
                end
                
                for i = 2:0 %length(mk)
                    included = 1;
                    for i2 = 1:length(mk{i})
                        found = 0;
                        for i3 = 1:i-1
                            bpms2 = mean([mk{i3}.lastbpm]);
                            norb = bpms2 ./ [mk{i3}.lvl];
                            nbpms = mean(norb) * [mk{i3}.lvl];
                            dist = min(abs(60/mean(mk{i}(i2).lastbpm) - 60./bpms2),...
                                       abs(60/mean(mk{i}(i2).lastbpm) - 60./nbpms));
                            if ~isempty(find(dist<option.lart));
                                found = 1;
                                %coord = [i2 i3];
                                break
                            end
                        end
                        if ~found
                            for i3 = 1:i-1
                                for i4 = 1:length(mk{i3})
                                    ma = max(mean(mk{i}(i2).lastbpm),...
                                             mean(mk{i3}(i4).lastbpm));
                                    mi = min(mean(mk{i}(i2).lastbpm),...
                                             mean(mk{i3}(i4).lastbpm));
                                    div = ma / mi;
                                    rdiv = round(div);
                                    if rdiv > 1 && ...
                                            ~isempty(find(~mod(rdiv,[2 3]))) && ...
                                            (mod(div,1) < option.lart2 || ...
                                             mod(div,1) > 1-option.lart2)
                                        if mean(mk{i}(i2).lastbpm) > ...
                                                mean(mk{i3}(i4).lastbpm)
                                            mk{i3}(end+1).lvl = ...
                                                mk{i}(i2).lvl * rdiv;
                                        else
                                            mk{i3}(end+1).lvl = ...
                                                mk{i}(i2).lvl / rdiv;
                                            if mk{i3}(end).lvl < 1
                                                for i5 = 1:length(mk{i3})
                                                    mk{i3}(i5).lvl = mk{i3}(i5).lvl * rdiv;
                                                end
                                            end
                                        end
                                        mk{i3}(end).lastbpm = mk{i}(i2).lastbpm;
                                        mk{i3}(end).bpms = mk{i}(i2).bpms;
                                        mk{i3}(end).timidx = mk{i}(i2).timidx;
                                        mk{i3}(end).score = mk{i}(i2).score;
                                        found = 1;
                                        bpms{i3}(:,end+1) = mk{i}(i2).lastbpm;
                                        %coord = [i2 i3];
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
                        mk(i) = [];
                        break
                    end
                end
            end
        end

        meters{j}{k} = mk;
    end
end

m = purgedata(p);
m = set(m,'Data',meters);
m = class(struct,'mirmetre',mirdata(m));
o = {m,p};


function best = findbest(m,timidx)
bestscore = -Inf;
best = [];
for i = 1:length(m)
    if m(i).timidx(end) == timidx && m(i).score(end) > bestscore
        best = i;
        bestscore = m(i).score(end);
    end
end


function bpm = getbpm(p,ptl)
if isa(p,'mirautocor') && not(get(p,'FreqDomain'))
    bpm = 60./ptl;
else
    bpm = ptl*60;
end