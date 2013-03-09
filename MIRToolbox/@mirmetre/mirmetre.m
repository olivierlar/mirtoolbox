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
        lart.default = .1;
    option.lart = lart;

        lart2.type = 'Integer';
        lart2.default = .2;
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
pv = get(p,'PeakVal');
for j = 1:length(pt)
    for k = 1:length(pt{j})
        ptk = pt{j}{k};
        for h = 1:size(ptk,3)
            mipv = +Inf;
            mapv = -Inf;
            for l = 1:length(pv{j}{k})
                if min(pv{j}{k}{l}) < mipv
                    mipv = min(pv{j}{k}{l});
                end
                if max(pv{j}{k}{l}) > mapv
                    mapv = max(pv{j}{k}{l});
                end
            end
            mk = {};
            globpm = [];
            for l = 1:size(ptk,2)       % For each successive frame
                %oldglobpm = globpm;
                ptl = getbpm(p,ptk{1,l,h}); % Peaks

                bpms = cell(1,length(mk));
                for i2 = 1:length(mk)
                    bpms{i2} = [mk{i2}.lastbpm];
                end
                
                foundk = zeros(1,length(mk));
                
                ampli = d{j}{k}(:,l,h);
                pos = ppp{j}{k}{l};
                                
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
                    
                    thri = (1-(pv{j}{k}{l}(i) - mipv)/(mapv - mipv))^2/10 ...
                           + .1;
                    
                    score = ampli(pos(i));
                    found = zeros(1,length(mk));              % Is peak in metrical hierarchies?
                    coord = [];             % Where is peak located
                    
                    dist = inf(1,length(mk));
                    indx = nan(1,length(mk));
                    dist2 = cell(1,length(mk));
                    for i2 = 1:length(mk)   % For each metrical hierarchy
                        if ~isempty(bpms{i2})
                            bpm2 = repmat(globpm(i2,end), [1 length(mk{i2})])...
                                   ./ [mk{i2}.lvl];
                            for i3 = 1:length(bpm2)
                                if 0 %mk{i2}(i3).timidx < l-5
                                    bpm2(i3) = NaN;
                                end
                            end
                            dist2{i2} = abs(60/ptli - 60./bpm2);
                            [dist(i2) indx(i2)] = min(dist2{i2});
                                
                            dist3 = abs(60/ptli - 60./[mk{i2}.lastbpm]);
                            [disti3 indx3] = min(dist3);
                            if disti3 < dist(i2) && indx(i2) ~= indx3 ...
                                    && ~mod(mk{i2}(indx3).lvl,1)
                                dist2{i2} = dist3;
                                dist(i2) = disti3;
                                indx(i2) = indx3;
                            end
                        end
                    end
                    [unused order] = sort(dist);
                    
                    for i2 = order
                        if foundk(i2)
                            thri2 = thri;
                        else
                            thri2 = min(thri,.05);
                        end
                        %locoord = [];
                        if dist(i2) < thri2
                            % Continuing an existing metrical level.
                            if mk{i2}(indx(i2)).timidx(end) ~= l
                                % Not already continued.
                                coord = [i2 indx(i2)];
                                if 1 %isempty(locoord) || ...
                                        %dist(i2) < dist2{i2}(locoord)
                                    % Level not identified to
                                    % one already detected
                                    mk{i2}(indx(i2)).timidx(end+1) = l;
                                    mk{i2}(indx(i2)).bpms(end+1) = ptl(i);
                                    mk{i2}(indx(i2)).lastbpm = ptli;
                                    mk{i2}(indx(i2)).score(end+1) = ...
                                        d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                    %locoord = indx(i2);

                                    if ~(foundk(i2))% && isempty(find(found,1))
                                        foundk(i2) = 1;
                                        mk{i2}(indx(i2)).major = 1;
                                    end
                                    if 0 %foundk(i2) == 1 && ...
                                            %~isempty(find(mk{i2}(indx(i2)).lvl ...
                                             %             == [.25 .5 1 2 3 4 5 6 8]))
                                        instnt = ptli * mk{i2}(indx(i2)).lvl;
                                        if d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h) > .1 ...
                                                && abs(60/instnt - 60/globpm(i2)) < .1
                                            globpm(i2) = instnt;
                                        end
                                        foundk(i2) = 2;
                                    end
                                end
                            end
                            found(i2) = 1;
                        end
                    end
                    
                    if 0 %~isempty(find(found))
                        continue
                    end
                    
                    incoherent = zeros(1,length(mk));
                    i2 = 1;
                    while i2 <= length(mk)
                        if found(i2)
                            i2 = i2+1;
                            continue
                        end
                        [unused ord] = sort(bpms{i2});
                        orbpms = bpms{i2}(ord);
                        fo = find(orbpms > ptli, 1);
                        if isempty(fo)
                            fo = size(orbpms,2)+1;
                        end

                        % Stored levels slower than candidate
                        slower = [];
                        i3 = fo-1;
                        err = Inf;
                        while i3 > 0
                            if l - mk{i2}(i3).timidx(end) > 10 || ...
                                mk{i2}(i3).score(end) < .1
                                i3 = i3-1;
                                continue
                            end
                                
                            bpm3 = globpm(i2,end) / mk{i2}(ord(i3)).lvl;
                            if ~foundk(i2)
                                div = [ptli ptli] ./ bpm3;
                            else
                                div = [ptli2; ptli1] ./ bpm3;
                            end
                            rdiv = round(ptli / bpm3);
                            if rdiv > 1 && rdiv < 7 %&& mk{i2}(ord(i3)).major
                                if floor(div(1)) ~= floor(div(2))
                                    newerr = 0;
                                else
                                    newerr = min(min(mod(div,1)),...
                                             min(1-mod(div,1)));
                                end
                                if 0 %~foundk(i2)
                                    thr = .02;
                                else
                                    thr = option.lart2;
                                end
                                if newerr < thr
                                    % Candidate level can be
                                    % integrated in this metrical
                                    % hierarchy
                                    if newerr < err
                                        if isempty(slower)
                                            slower.lvl = mk{i2}(ord(i3)).lvl / rdiv;
                                            slower.bpm = orbpms(i3);
                                            slower.score = mk{i2}(ord(i3)).score(end);
                                            slower.rdiv = rdiv;
                                            rptli1 = orbpms(i3) * (rdiv - .4);
                                            if ptli1 < rptli1
                                                ptli1 = rptli1;
                                            end
                                            rptli2 = orbpms(i3) * (rdiv + .4);
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
                        while i3 <= length(orbpms)
                            if l - mk{i2}(i3).timidx(end) > 10
                                i3 = i3+1;
                                continue
                            end
                            
                            bpm3 = globpm(i2,end) / mk{i2}(ord(i3)).lvl;
                            if ~foundk(i2)
                                div = bpm3 ./ [ptli ptli];
                            else
                                div = bpm3 ./ [ptli2;ptli1];
                            end
                            rdiv = round(bpm3 / ptli);
                            if rdiv > 1 % && rdiv < 7 %&& mk{i2}(ord(i3)).major
                                if floor(div(1)) < floor(div(2))
                                    newerr = 0;
                                else
                                    newerr = min(min(mod(div,1)),...
                                             min(1-mod(div,1)));
                                end
                                if ~foundk(i2)
                                    thr = .01;
                                else
                                    thr = option.lart2;
                                end
                                if newerr < thr
                                    % Candidate level can be
                                    % integrated in this metrical
                                    % hierarchy
                                    if newerr < err
                                        if isempty(faster)
                                            faster.lvl = mk{i2}(ord(i3)).lvl * rdiv;
                                            faster.bpm = orbpms(i3);
                                            faster.score = mk{i2}(ord(i3)).score(end);
                                            faster.rdiv = rdiv;
                                            rptli1 = orbpms(i3) / (rdiv + .4);
                                            if ptli1 < rptli1
                                                ptli1 = rptli1;
                                            end
                                            rptli2 = orbpms(i3) / (rdiv - .4);
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
                                                    
                        found(i2) = 1;
                        l0 = find(lvl == [mk{i2}.lvl]);
                        if isempty(l0)
                            % New metrical level
                            mk{i2}(end+1).lvl = lvl;
                            mk{i2}(end).lastbpm = ptli;
                            mk{i2}(end).bpms = ptl(i);
                            mk{i2}(end).timidx = l;
                            mk{i2}(end).score = ampli(pos(i));
                            mk{i2}(end).major = 0;
                            coord = [i2 length(mk{i2})];
                            bpms{i2}(end+1) = ptli;
                        else
                            dist = abs(mean([mk{i2}(l0).lastbpm]) - ptl(i));
                            [unused md] = min(dist);
                            if mk{i2}(l0(md)).timidx(end) < l
                                mk{i2}(l0(md)).lastbpm = ptli;
                                mk{i2}(l0(md)).bpms(end+1) = ptl(i);
                                mk{i2}(l0(md)).score(end+1) = ampli(pos(i));
                                mk{i2}(l0(md)).timidx(end+1) = l;
                            elseif score > mk{i2}(l0(md)).score
                                mk{i2}(l0(md)).lastbpm = ptli;
                                mk{i2}(l0(md)).bpms(end) = ptl(i);
                                mk{i2}(l0(md)).score(end) = ampli(pos(i));
                            end
                            coord = [i2 l0(md)];
                        end
                        %locoord = coord(2);
                        if ~(foundk(i2))
                            foundk(i2) = 1;
                            mk{i2}(coord(2)).major = 1;
                        end
                        if 0 %foundk(i2) == 1 && ...
                             %   ~isempty(find(mk{i2}(coord(2)).lvl ...
                            %                  == [.25 .5 1 2 3 4 5 6 8]))
                            instnt = ptli * lvl;
                            if abs(60/instnt - 60/globpm(i2)) < .01
                                globpm(i2) = instnt;
                            end
                            foundk(i2) = 2;
                        end
                        
                        i2 = i2 + 1;
                    end
                                        
                    if isempty(mk) && isempty(find(found)) && isempty(find(incoherent)) ...
                            && d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h) > .15
                        % New metrical hierarchy
                        mk{end+1}.lvl = 1;
                        mk{end}.lastbpm = ptli;
                        mk{end}.bpms = ptl(i);
                        mk{end}.timidx = l;
                        mk{end}.score = ...
                            d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                        mk{end}.globpms = [];
                        mk{end}.major = 1;
                        %found(end+1) = 1;
                        bpms{end+1} = ptli;
                        foundk(end+1) = 1;
                        globpm(end+1,:) = NaN(1,l);
                        globpm(end,l) = ptli;
                        %coord = [length(bpms),1];
                    end
                end
                
                for i = 1:length(mk)
                    if l == 1 || isnan(globpm(i,l-1))
                        glo = 0;
                        sco = 0;
                        for i2 = 1:length(mk{i})
                            if mk{i}(i2).timidx(end) == l
                                sco2 = mk{i}(i2).score(end);
                                ind = mk{i}(i2).bpms(end) * mk{i}(i2).lvl;
                                glo = glo + ind * sco2;
                                sco = sco + sco2;
                            end
                        end
                        globpm(i,l) = glo / sco;
                    else
                        glog = log2(globpm(i,l-1));
                        glodif = 0;
                        sco = 0;
                        mindif = Inf;
                        for i2 = 1:length(mk{i})
                            if mk{i}(i2).timidx(end) == l
                                sco2 = mk{i}(i2).score(end);
                                dif = glog - log2(mk{i}(i2).bpms(end) ...
                                                  * mk{i}(i2).lvl);
                                glodif = glodif + dif * sco2;
                                sco = sco + sco2;
                                if abs(dif) < abs(mindif)
                                    mindif = dif;
                                end
                            end
                        end
                        if glodif
                            glodif = glodif / sco;
                        end
                        if abs(glodif) > abs(mindif)
                            if glodif * mindif < 0
                                glodif = 0;
                            else
                                glodif = mindif;
                            end
                        end
                        globpm(i,l) = globpm(i,l-1) / 2^glodif;
                    end
                    for i2 = 1:length(mk{i})
                        if mk{i}(i2).timidx(end) == l
                            mk{i}(i2).globpms(end+1) = globpm(i,l) ...
                                                        / mk{i}(i2).lvl;
                        end
                    end
                end
                
                for i = 2:length(mk)
                    included = 1;
                    for i2 = 1:length(mk{i})
                        found = 0;
                        for i3 = 1:i-1
                            nbpms1 = globpm(i,l)/ mk{i}(i2).lvl;
                            nbpms2 = repmat(globpm(i3,l),[1,size(mk{i3},2)])...
                                                ./ [mk{i3}.lvl];
                            dist = abs(60/nbpms1 - 60./nbpms2);
                            i4 = find(dist<.1,1);
                            if ~isempty(i4);
                                found = 1;
                                break
                                ratio = mk{i3}(i4).lvl / mk{i}(i2).lvl;
                                for i4 = 1:length(mk{i})
                                    lvl = mk{i}(i4).lvl * ratio;
                                    if isempty(find(lvl == [mk{i3}.lvl],1))
                                        mk{i3}(end+1).lvl = lvl;
                                        mk{i3}(end).lastbpm = mk{i}(i4).lastbpm;
                                        mk{i3}(end).bpms = mk{i}(i4).bpms(end);
                                        mk{i3}(end).globpms = globpm(i,l) ...
                                                    / lvl;
                                        mk{i3}(end).timidx = mk{i}(i2).timidx(end);
                                        mk{i3}(end).score = mk{i}(i2).score(end);
                                        found = 1;
                                        bpms{i3}(end+1) = mk{i}(i4).lastbpm;
                                    end
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
                        globpm(i,:) = [];
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


function bpm = getbpm(p,ptl)
if isa(p,'mirautocor') && not(get(p,'FreqDomain'))
    bpm = 60./ptl;
else
    bpm = ptl*60;
end