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
        frame.keydefault = [3 .1];
    option.frame = frame;
    
        fea.type = 'String';
        fea.choice = {'Envelope','DiffEnvelope','SpectralFlux',...
                      'Pitch','Novelty'};
        fea.default = 'Novelty';
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
y = mirautocor(x,'Min',60/option.ma,'Max',60/option.mi,...
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
            oldmeters = {};
            for l = 1:size(ptk,2)       % For each successive frame
                ptl = getbpm(p,ptk{1,l,h}); % Peaks
                cntr = zeros(1,length(mk));
                
                bpms = cell(1,length(mk));
                for i2 = 1:length(mk)
                    bpms{i2} = [mk{i2}.lastbpm];
                end
                
                for i = 1:length(ptl)   % For each peak
                    found = 0;          % Is peak located in metrical hierarchies?
                    coord = [];         % Where is peak located
                    for i2 = 1:length(mk)   
                                        % For each metrical hierarchy
                        if ~isempty(bpms{i2})
                            locoord = [];
                            dist = abs(60/ptl(i) - 60./bpms{i2});

                            i3 = 1;
                            while i3 <= length(dist)
                                if dist(i3) < option.lart
                                    % Continuing an existing metrical
                                    % level.
                                    if mk{i2}(i3).timidx(end) == l
                                        % Already continued.
                                        odist = abs(60./mk{i2}(i3).bpms(end) ...
                                                    - 60./bpms{i2}(i3));
                                        if dist(i3) < odist
                                            % New candidate is better.
                                            mk{i2}(i3).bpms(end) = ptl(i);
                                            mk{i2}(i3).score(end) = ...
                                                d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                        end

                                        locoord = i3;
                                    else
                                        if bpms{i2}(i3) > ptl(i)
                                            div = bpms{i2}(i3) / ptl(i);
                                        else
                                            div = ptl(i) / bpms{i2}(i3);
                                        end
                                        if div < 1+option.lart2
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
                                                mk{i2}(i3).lastbpm = ptl(i);
                                                mk{i2}(i3).score(end+1) = ...
                                                    d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);

                                                locoord = i3;
                                            end
                                        end
                                    end
                                end
                                i3 = i3 + 1;
                            end
                        end
                    end

                    if found
                        if 0 %~isempty(locoord)
                            cntr(i2) = cntr(i2)+1;
                            if cntr(i2) ~= locoord
                                indx = 1:length(meters{i2});
                                indx(indx == locoord) = [];
                                indx = [indx(1:cntr(i2)-1) ...
                                        locoord ...
                                        indx(cntr(i2):end)];
                                meters{i2} = meters{i2}(indx);
                            end
                        end
                    else
                        % Candidate level not belonging to current
                        % metrical levels
                        i2 = 1;
                        while i2 <= length(mk)
                            locoord = [];
                            i3 = 1;
                            while i3 <= length(bpms{i2})
                                if ptl(i) > bpms{i2}(i3)
                                    % Candidate faster than stored
                                    % level
                                    div = ptl(i) / bpms{i2}(i3);
                                    rdiv = round(div);
                                    if rdiv > 1 && ...
                                            ...~isempty(find(~mod(rdiv,[2 3 5]))) && ...
                                            (mod(div,1) < option.lart2 || ...
                                             mod(div,1) > 1-option.lart2)
                                        % Candidate level can be
                                        % integrated in this metrical
                                        % hierarchy

                                        if ~isempty(locoord)
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

                                        elseif found
                                            % Candidate level also
                                            % integrated in other
                                            % metrical hierarchy. Both
                                            % hierarchies are fused.
                                            lvl = mk{i2}(i3).lvl / rdiv;
                                            if lvl<1
                                                for i4 = 1:length(mk{i2})
                                                    mk{i2}(i4).lvl = mk{i2}(i4).lvl * rdiv;
                                                end
                                                lvl = mk{i2}(i3).lvl / rdiv;
                                            end
                                            chcur = 0;
                                            if mk{coord(1)}(coord(2)).lvl > ...
                                                    mk{i2}(i3).lvl / rdiv
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
                                                lvl1 = lvl;
                                                    % Level in slower
                                                lvl2 = meter2(coord(2)).lvl;
                                                    % Level in faster
                                            else
                                                % Other hierarchy is
                                                % slower than current
                                                meter1 = mk{coord(1)};
                                                meter2 = mk{i2};
                                                bpms{coord(1)} = bpms{i2};
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
                                                    meter2(i4).lvl = meter2(i4).lvl * mult2;
                                                end
                                            end
                                            for i4 = 1:length(meter1)
                                                % Fusing one hierarchy
                                                % into the other..
                                                f = find(meter1(i4).lvl * mult1 ...
                                                             == [meter2.lvl]);
                                                if isempty(f)
                                                    meter2(end+1) = meter1(i4);
                                                    meter2(end).lvl = meter2(end).lvl * mult1;
                                                    f = length(meters);
                                                    bpms{coord(1)}(end+1) = meter1(i4).lastbpm;
                                                end
                                                coord(2) = f;
                                            end
                                            mk{coord(1)} = meter2;
                                            mk(i2) = [];
                                            bpms(i2) = [];
                                            i2 = i2 - 1;
                                            break
                                        else
                                            found = 1;
                                            lvl = mk{i2}(i3).lvl ...
                                                  / round(div);
                                            l0 = find(lvl == ...
                                                      [mk{i2}.lvl]);
                                            score = d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                            if isempty(l0)
                                                % New metrical level
                                                mk{i2}(end+1).lvl = lvl;
                                                mk{i2}(end).lastbpm = ptl(i);
                                                mk{i2}(end).bpms = ptl(i);
                                                mk{i2}(end).timidx = l;
                                                mk{i2}(end).score = score;
                                                coord = [i2 length(mk{i2})];
                                                bpms{i2}(end+1) = ptl(i);
                                                if lvl<1
                                                    for i4 = 1:length(mk{i2})
                                                        mk{i2}(i4).lvl = mk{i2}(i4).lvl * round(div);
                                                    end
                                                end
                                            else
                                                dist = abs([mk{i2}(l0).lastbpm] - ptl(i));
                                                [unused md] = min(dist);
                                                if mk{i2}(l0(md)).timidx(end) < l
                                                    mk{i2}(l0(md)).lastbpm = ptl(i);
                                                    mk{i2}(l0(md)).bpms(end+1) = ptl(i);
                                                    mk{i2}(l0(md)).score(end+1) = score;
                                                    mk{i2}(l0(md)).timidx(end+1) = l;
                                                elseif score > mk{i2}(l0(md)).score
                                                    mk{i2}(l0(md)).lastbpm = ptl(i);
                                                    mk{i2}(l0(md)).bpms(end) = ptl(i);
                                                    mk{i2}(l0(md)).score(end) = score;
                                                end
                                                coord = [i2 l0(md)];
                                            end
                                        end
                                        locoord = coord(2);
                                    end
                                else
                                    % Candidate slower than stored
                                    % level
                                    div = bpms{i2}(i3) / ptl(i);
                                    rdiv = round(div);
                                    if rdiv > 1 && ...
                                            ...~isempty(find(~mod(rdiv,[2 3 5]))) && ...
                                            (mod(div,1) < option.lart2 || ...
                                             mod(div,1) > 1-option.lart2)
                                        % Candidate level can be
                                        % integrated in this metrical
                                        % hierarchy

                                        if ~isempty(locoord)
                                            % Level identified to one
                                            % already detected

                                            %meters{i2}(i3) = [];
                                            %bpms(i3) = [];
                                            %i3 = i3 - 1;
                                            %if res(1) == i2 && res(2)>i3
                                            %    res(2) = res(2)-1;
                                            %    res3 = res;
                                            %end
                                        else
                                            score = d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                                            lvl = mk{i2}(i3).lvl ...
                                                    * round(div);
                                            if 1 %lvl == round(lvl)
                                                l0 = find(lvl == ...
                                                          [mk{i2}.lvl]);
                                                if isempty(l0)
                                                    % New metrical level
                                                    mk{i2}(end+1).lvl = lvl;
                                                    mk{i2}(end).lastbpm = ptl(i);
                                                    mk{i2}(end).bpms = ptl(i);
                                                    mk{i2}(end).timidx = l;
                                                    mk{i2}(end).score = score;
                                                    found = 1;
                                                    coord = [i2 length(mk{i2})];
                                                    bpms{i2}(end+1) = ptl(i);
                                                    %comet(i2) = 1;
                                                else
                                                    dist = abs([mk{i2}(l0).lastbpm] - ptl(i));
                                                    [unused md] = min(dist);
                                                    if mk{i2}(l0(md)).timidx(end) < l
                                                        mk{i2}(l0(md)).lastbpm = ptl(i);
                                                        mk{i2}(l0(md)).bpms(end+1) = ptl(i);
                                                        mk{i2}(l0(md)).score(end+1) = score;
                                                        mk{i2}(l0(md)).timidx(end+1) = l;
                                                        found = 1;
                                                        coord = [i2 l0(md)];
                                                    elseif score > mk{i2}(l0(md)).score
                                                        mk{i2}(l0(md)).lastbpm = ptl(i);
                                                        mk{i2}(l0(md)).bpms(end) = ptl(i);
                                                        mk{i2}(l0(md)).score(end) = score;
                                                        found = 1;
                                                        coord = [i2 0];
                                                    end
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
                                                coord = [i2 length(meters{i2})];
                                                comet(i2) = 1;
                                                bpms{i2}(end+1) = ptl(i);
                                            end
                                        end
                                        if ~isempty(coord)
                                            locoord = coord(2);
                                        end
                                    end
                                end
                                i3 = i3 + 1;
                            end
                            i2 = i2 + 1;
                        end
                        %for i2 = 1:length(meters)
                        %    if ~isempty(meters{i2})
                        %        [unused ord] = sort([meters{i2}.lastbpm],'descend');
                        %        meters{i2} = meters{i2}(ord);
                        %    end
                        %end

                        if found
                            if 0 %coord(2)
                                cntr(coord(1)) = cntr(coord(1))+1;
                                if cntr(coord(1)) ~= coord(2)
                                    indx = 1:length(meters{coord(1)});
                                    indx(indx == coord(2)) = [];
                                    indx = [indx(1:cntr(coord(1))-1) ...
                                            coord(2) ...
                                            indx(cntr(coord(1)):end)];
                                    meters{coord(1)} = meters{coord(1)}(indx);
                                end
                            end
                        else
                            % New metrical hierarchy
                            mk{end+1}.lvl = 1;
                            mk{end}.lastbpm = ptl(i);
                            mk{end}.bpms = ptl(i);
                            mk{end}.timidx = l;
                            mk{end}.score = ...
                                d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                            cntr(end+1) = 0;
                            %found = 1;
                            bpms{end+1} = ptl(i);
                            %coord = [length(bpms),1];
                        end
                    end
                end

                for i = 2:length(mk)
                    if length(oldmeters)<i || isempty(oldmeters{i})
                        included = 1;
                        for i2 = 1:length(mk{i})
                            found = 0;
                            for i3 = 1:i-1
                                bpms2 = [mk{i3}.lastbpm];
                                norb = bpms2 ./ [mk{i3}.lvl];
                                nbpms = mean(norb) * [mk{i3}.lvl];
                                dist = min(abs(60/mk{i}(i2).lastbpm - 60./bpms2),...
                                           abs(60/mk{i}(i2).lastbpm - 60./nbpms));
                                if ~isempty(find(dist<option.lart));
                                    found = 1;
                                    %coord = [i2 i3];
                                    break
                                end
                            end
                            if ~found
                                for i3 = 1:i-1
                                    for i4 = 1:length(mk{i3})
                                        ma = max(mk{i}(i2).lastbpm,...
                                                 mk{i3}(i4).lastbpm);
                                        mi = min(mk{i}(i2).lastbpm,...
                                                 mk{i3}(i4).lastbpm);
                                        div = ma / mi;
                                        rdiv = round(div);
                                        if rdiv > 1 && ...
                                                ~isempty(find(~mod(rdiv,[2 3]))) && ...
                                                (mod(div,1) < option.lart2 || ...
                                                 mod(div,1) > 1-option.lart2)
                                            if mk{i}(i2).lastbpm > ...
                                                    mk{i3}(i4).lastbpm
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
                                            bpms{i3}(end+1) = mk{i}(i2).lastbpm;
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
                            bpms2(i) = [];
                            if length(oldmeters) >= i
                                oldmeters(i) = [];
                            end
                            break
                        end
                    end
                end

                for i = 1:length(mk)
                    i2 = 1;
                    while i2 <= length(mk{i})
                        if mk{i}(i2).timidx(end) ~= l || ...
                                l == size(ptk,2)
                            if length(oldmeters) < i
                                oldmeters{i} = mk{i}(i2);
                            else
                                oldmeters{i}(end+1) = mk{i}(i2);
                            end
                            mk{i}(i2) = [];
                            bpms{i}(i2) = [];
                        else
                            mk{i}(i2).lastbpm = mk{i}(i2).bpms(end);
                            i2 = i2+1;
                        end
                    end
                end
                if ~mod(l,100)
                    l
                end
            end
        end

        mk = oldmeters;
        meters{j}{k} = mk;
    end
end

m.meters = meters;
m = class(m,'mirmetre');
o = {m,p};


function bpm = getbpm(p,ptl)
if isa(p,'mirautocor') && not(get(p,'FreqDomain'))
    bpm = 60./ptl;
else
    bpm = ptl*60;
end