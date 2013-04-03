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
    
        envelope.key = 'Envelope';
        envelope.type = 'Boolean';
        envelope.default = 0;
    option.envelope = envelope;
    
        novelty.key = 'Novelty';
        novelty.type = 'Boolean';
        novelty.default = 0;
    option.novelty = novelty;
    
        mix.key = 'Mix';
        mix.type = 'Boolean';
        mix.default = 0;
    option.mix = mix;
    
        sgate.key = 'SmoothGate';
        sgate.type = 'Boolean';
        sgate.default = 0;
    option.sgate = sgate;
    
%% option related to mirautocor:                
    
        nw.key = 'NormalWindow';
        nw.default = 0;
    option.nw = nw;

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

if isamir(x,'mirautocor')
    y = mirautocor(x,'Min',60/option.ma,'Max',60/option.mi * 2,...
          'NormalWindow',option.nw);
    if ischar(option.sum)
        y = mirsum(y);
    end
elseif isamir(x,'mirenvelope')
    x = mironsets(x,'Frame',option.frame.length.val,...
                            option.frame.length.unit,...
                            option.frame.hop.val,...
                            option.frame.hop.unit);
    y = mirautocor(x,'Min',60/option.ma,'Max',60/option.mi * 2,...
                     'NormalWindow',option.nw);
elseif option.mix
    o1 = mironsets(x,'Diff','Detect',0,...
                    'Frame',option.frame.length.val,...
                            option.frame.length.unit,...
                            option.frame.hop.val,...
                            option.frame.hop.unit);
    ac1 = mirautocor(o1,'Min',60/option.ma,'Max',60/option.mi * 2,...
          'NormalWindow',option.nw);

    o2 = mironsets(x,'Diff','Novelty','Detect',0,...
                    'Frame',option.frame.length.val,...
                            option.frame.length.unit,...
                            option.frame.hop.val,...
                            option.frame.hop.unit);
    ac2 = mirautocor(o2,'Min',60/option.ma,'Max',60/option.mi * 2,...
          'NormalWindow',option.nw);
    y = max(ac1,ac2);

elseif option.novelty
    o = mironsets(x,'Diff','Novelty','Detect',0,...
                    'Frame',option.frame.length.val,...
                            option.frame.length.unit,...
                            option.frame.hop.val,...
                            option.frame.hop.unit);
    y = mirautocor(o,'Min',60/option.ma,'Max',60/option.mi * 2,...
          'NormalWindow',option.nw);
elseif option.sgate
    o = mironsets(x,'SmoothGate','Detect',0,...
                    'Frame',option.frame.length.val,...
                            option.frame.length.unit,...
                            option.frame.hop.val,...
                            option.frame.hop.unit);
    y = mirautocor(o,'Min',60/option.ma,'Max',60/option.mi * 2,...
          'NormalWindow',option.nw);
else
    o = mironsets(x,'Diff','Detect',0,...
                    'Frame',option.frame.length.val,...
                            option.frame.length.unit,...
                            option.frame.hop.val,...
                            option.frame.hop.unit);
    y = mirautocor(o,'Min',60/option.ma,'Max',60/option.mi * 2,...
          'NormalWindow',option.nw);

end

y = mirpeaks(y,'Total',Inf,...
               'Threshold',option.thr,'Contrast',option.cthr,...
               'NoBegin','NoEnd',...
               'Normalize','Local','Order','Amplitude');
type = {'mirmetre',mirtype(y)};
    

function m = main(p,option,postoption)
if iscell(p)
    p = p{1};
end
if isamir(p,'mirscalar')
    m = modif(m,postoption);
    m = {t};
    return
end
pt = get(p,'PeakPrecisePos');
meters = cell(1,length(pt));
globpms = cell(1,length(pt));
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
                    for i2 = 1:length(mk)   % For each metrical hierarchy
                        if ~mk{i2}(1).active || isempty(bpms{i2})
                            continue
                        end
                        
                        %if l > 1 && ~isnan(globpm(i2,l-1))
                        %    globpmi2 = globpm(i2,l-1);
                        %else
                        %    globpmi2 = globpm(i2,l);
                        %end
                        globpmi2 = globpm(i2,end);
                        
                        bpm2 = repmat(globpmi2, [1 length(mk{i2})])...
                               ./ [mk{i2}.lvl];
                        dist2 = abs(60/ptli - 60./bpm2);
                        if 1 %~foundk(i2)
                            nf = [];
                            for i3 = 1:length(mk{i2})
                                if ~isempty(mk{i2}(i3).function)
                                    nf(end+1) = i3;
                                end
                            end
                            lvls = sort([mk{i2}(nf).lvl]);
                            for i3 = 1:length(mk{i2})
                                lvl3 = mk{i2}(i3).lvl;
                                if isempty(mk{i2}(i3).function)
                                    l1 = find(lvls(end:-1:1) < mk{i2}(i3).lvl,1);
                                    l2 = find(lvls > mk{i2}(i3).lvl,1);
                                    if (~isempty(l1) && ...
                                            mod(lvl3,lvls(l1))) || ...
                                            (~isempty(l2) && ...
                                             mod(lvls(l2),lvl3)) 
                                        dist2(i3) = NaN;
                                    end
                                end
                            end
                        end
                        [disti2 indx2] = min(dist2);
                        
                        if 0
                            dist2b = dist2;
                            for i3 = 1:length(mk{i2})
                                if isempty(mk{i2}(i3).function)
                                    dist2b(i3) = NaN;
                                end
                            end
                            [disti2b indx2b] = min(dist2b);
                            if foundk(i2)
                                thri2 = thri;
                            else
                                thri2 = min(thri,.1);
                            end
                            if indx2b ~= indx2 && disti2b < thri2
                                disti2 = disti2b;
                                indx2 = indx2b;
                            end
                        end

                        dist3 = NaN(1,length(mk{i2}));
                        for i3 = 1:length(mk{i2})
                            t3 = find(mk{i2}(i3).timidx == l-1);
                            if ~isempty(t3)
                                dist3(i3) = abs(60/ptli - ...
                                                60./mk{i2}(i3).bpms(t3));
                            end
                        end
                        [disti3 indx3] = min(dist3);
                        if disti3 < disti2
                            dist(i2) = disti3;
                            indx(i2) = indx3;
                        else
                            dist(i2) = disti2;
                            indx(i2) = indx2;
                            if abs(log2(ptli / bpm2(indx2))) > .3
                                dist(i2) = Inf;
                                indx(i2) = 0;
                                continue
                            end
                        end
                    end
                    
                    [unused order] = sort(dist);
                    for i2 = order
                        if foundk(i2)
                            thri2 = thri;
                        else
                            thri2 = min(thri,.1);
                        end
                        
                        if dist(i2) > thri2
                            continue
                        end
                        
                        % Continuing an existing metrical level.
                        if mk{i2}(indx(i2)).timidx(end) ~= l
                            coord = [i2 indx(i2)];
                            % Level not identified to
                            % one already detected
                            mk{i2}(indx(i2)).timidx(end+1) = l;
                            mk{i2}(indx(i2)).bpms(end+1) = ptl(i);
                            mk{i2}(indx(i2)).lastbpm = ptli;
                            mk{i2}(indx(i2)).score(end+1) = ...
                                d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);

                            if ~foundk(i2)% && isempty(find(found,1))
                                if isempty(find(foundk,1)) && ...
                                        isempty(mk{i2}(indx(i2)).function) && ...
                                         mk{i2}(indx(i2)).score(end) > .3
                                    if mk{i2}(indx(i2)).reldiv > 0
                                        ref = mk{i2}(indx(i2)).lvl / mk{i2}(indx(i2)).reldiv;
                                    else
                                        ref = mk{i2}(indx(i2)).lvl * -mk{i2}(indx(i2)).reldiv;
                                    end
                                    mk{i2}(indx(i2)).function = ...
                                        [reldiv; ref];
                                end
                                foundk(i2) = 1;
                                if ~isempty(mk{i2}(indx(i2)).function)
                                    globpm(i2,l) = ptli * mk{i2}(indx(i2)).lvl;
                                    for i3 = 1:size(globpm,1)
                                        if globpm(i3,l) == 0
                                            globpm(i3,l) = globpm(i3,l-1);
                                        end
                                    end
                                end
                            end
                        elseif abs(mk{i2}(indx(i2)).bpms(end) - ...
                                   globpm(i2,end) / mk{i2}(indx(i2)).lvl) > ...
                               abs(ptl(i) - ...
                                   globpm(i2,end) / mk{i2}(indx(i2)).lvl)
                            coord = [i2 indx(i2)];
                            % Level already identified to
                            % one already detected
                            mk{i2}(indx(i2)).bpms(end) = ptl(i);
                            mk{i2}(indx(i2)).lastbpm = ptli;
                            mk{i2}(indx(i2)).score(end) = ...
                                d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                        end
                        found(i2) = 1;
                    end
                    
                    i2 = 1;
                    while i2 <= length(mk)
                        if found(i2)
                            i2 = i2+1;
                            continue
                        end
                        if ~mk{i2}(1).active
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
                            if l - mk{i2}(ord(i3)).timidx(end) > 10 || ...
                                mk{i2}(ord(i3)).score(end) < .1
                                i3 = i3-1;
                                continue
                            end
                            
                            if mk{i2}(ord(i3)).reldiv > 0
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
                            if rdiv == 0 || rdiv >= 7
                                i3 = i3-1;
                                continue
                            end
                            lvl = mk{i2}(ord(i3)).lvl / rdiv;
                            
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
                            
                            if newerr > thr
                                i3 = i3-1;
                                continue
                            end
                            
                            if 0 && ~isempty(find([mk{i2}.lvl] > lvl & ...
                                             [mk{i2}.lvl] < mk{i2}(ord(i3)).lvl))
                                i3 = i3-1;
                                continue
                            end
                            
                            % Candidate level can be
                            % integrated in this metrical
                            % hierarchy
                            if newerr < err
                                if isempty(slower)
                                    slower.ref = ord(i3);
                                    slower.lvl = lvl;
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
                                    %ptli = mean([ptli1,ptli2]);
                                    err = newerr;
                                    %break
                                elseif mk{i2}(ord(i3)).lvl / rdiv ...
                                        ~= slower.lvl
                                    slower.ref = ord(i3);
                                    slower.lvl = lvl;
                                    slower.rdiv = rdiv;
                                    %slower = [];
                                    %break
                                end
                            end
                            
                            i3 = i3 - 1;
                        end

                        % Stored levels faster than candidate
                        faster = [];
                        i3 = fo;
                        err = Inf;
                        while i3 <= length(orbpms)
                            if l - mk{i2}(ord(i3)).timidx(end) > 10
                                i3 = i3+1;
                                continue
                            end
                            
                            %if ...~foundk(i2) &&
                            %        isempty(mk{i2}(ord(i3)).function)
                            %    i3 = i3+1;
                            %    continue
                            %end
                            
                            bpm3 = globpm(i2,end) / mk{i2}(ord(i3)).lvl;
                            if ~foundk(i2)
                                div = bpm3 ./ [ptli ptli];
                            else
                                div = bpm3 ./ [ptli2;ptli1];
                            end
                            rdiv = round(bpm3 / ptli);
                            if rdiv > 1
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
                                            faster.ref = ord(i3);
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
                                            faster.ref = ord(i3);
                                            faster.lvl = mk{i2}(ord(i3)).lvl * rdiv;
                                            faster.rdiv = rdiv;
                                            %faster = [];
                                            %break
                                        end
                                    end
                                end
                            end
                            i3 = i3 + 1;
                        end

                        if isempty(slower) && isempty(faster)
                            i2 = i2 + 1;
                            continue
                        elseif isempty(slower)
                            lvl = faster.lvl;
                            rdiv = faster.rdiv;
                            reldiv = rdiv;
                            ref = faster.ref;
                        elseif isempty(faster)
                            lvl = slower.lvl;
                            rdiv = slower.rdiv;
                            reldiv = -rdiv;
                            ref = slower.ref;
                        elseif slower.score < faster.score
                            lvl = faster.lvl;
                            rdiv = faster.rdiv;
                            reldiv = rdiv;
                            ref = faster.ref;
                        else
                            lvl = slower.lvl;
                            rdiv = slower.rdiv;
                            reldiv = -rdiv;
                            ref = slower.ref;
                        end
                                                    
                        found(i2) = 1;
                        l0 = find(lvl == [mk{i2}.lvl]);
                        if isempty(l0)
                            % New metrical level
                            mk{i2}(end+1).lvl = lvl;
                            if isempty(mk{i2}(ref).function)
                                mk{i2}(end).function = [];
                            else
                                same = find(mk{i2}(ref).function(1,:)...
                                            *reldiv < 0);
                                saillant = 0; %isempty(same);
                                for i3 = 1:length(same)
                                    refdiv = mk{i2}(ref)...
                                                .function(1,same(i3));
                                    otherlvl = mk{i2}(ref)...
                                                .function(2,same(i3));
                                    other = find([mk{i2}.lvl] ...
                                                 == otherlvl);
                                    if abs(reldiv) < abs(refdiv)
                                        intradiv = abs(refdiv/reldiv);
                                        if ~isempty(mk{i2}(other).function)
                                            otherfunction = ...
                                                find(mk{i2}(other).function(1,:)...
                                                     == -refdiv);
                                            if round(intradiv) == intradiv
                                                mk{i2}(ref)...
                                                    .function(:,same(i3)) = ...
                                                        [-reldiv; lvl];
                                                mk{i2}(other).function(:,otherfunction)...
                                                    = [intradiv; lvl];
                                            end
                                        end
                                        mk{i2}(end).function = ...
                                            [reldiv,-intradiv; ...
                                             mk{i2}(ref).lvl,...
                                             mk{i2}(other).lvl];
                                        saillant = 2;
                                        break
                                    elseif mk{i2}(other).timidx(end)...
                                                ~= l || ...
                                            mk{i2}(other).score(end)...
                                                < ampli(pos(i))
                                        saillant = 1;
                                    end
                                end
                                if saillant == 1
                                    mk{i2}(end).function = ...
                                        [reldiv; mk{i2}(ref).lvl];
                                    mk{i2}(ref).function(:,end+1) = ...
                                        [-reldiv; lvl];
                                elseif saillant == 0
                                    mk{i2}(end).function = [];
                                end
                            end
                            mk{i2}(end).lastbpm = ptli;
                            mk{i2}(end).bpms = ptl(i);
                            mk{i2}(end).timidx = l;
                            mk{i2}(end).score = ampli(pos(i));
                            mk{i2}(end).ref = mk{i2}(ref).lvl;
                            mk{i2}(end).reldiv = reldiv;
                            
                            if reldiv < 0 && ~isempty(mk{i2}(ref).element)
                                mk{i2}(ref).element = [];
                                mk{i2}(end).element = 1;
                            end
                            
                            for i3 = 1:length(mk{i2})
                                if ptli > mk{i2}(i3).lastbpm
                                    div = ptli / mk{i2}(i3).lastbpm;
                                    if round(div) > 1 && ...
                                            min(mod(div,1),1-mod(div,1)) < .2
                                        ptli3 = mk{i2}(i3).lastbpm;
                                        nlvl = lvl * round(div);
                                        olvl = mk{i2}(i3).lvl;
                                        if nlvl ~= olvl && ...
                                                abs(60/ptli3 - ...
                                                    60/(globpm(i2,l)/nlvl)) < ...
                                                abs(60/ptli3 - ...
                                                    60/(globpm(i2,l)/olvl))
                                            mk{i2}(i3).lvl = nlvl;
                                            mk{i2}(i3).function = [];
                                            mk{i2}(i3).ref = lvl;
                                            mk{i2}(i3).reldiv = round(div);
                                        end
                                    end
                                end
                            end
                            
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
                        if ~(foundk(i2))
                            foundk(i2) = 1;
                            globpm(i2,l) = ptli * mk{i2}(coord(2)).lvl;
                            for i3 = 1:size(globpm,1)-1
                                if globpm(i3,l) == 0
                                    globpm(i3,l) = globpm(i3,l-1);
                                end
                            end
                        end
                        
                        i2 = i2 + 1;
                    end
                                        
                    if isempty(find(found))...
                            && d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h) > .15
                        % New metrical hierarchy
                        mk{end+1}.lvl = 1;
                        mk{end}.function = [0;1];
                        mk{end}.lastbpm = ptli;
                        mk{end}.bpms = ptl(i);
                        mk{end}.timidx = l;
                        mk{end}.score = ...
                            d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);
                        mk{end}.globpms = [];
                        mk{end}.locked = 0;
                        mk{end}.active = 1;
                        mk{end}.ref = [];
                        mk{end}.element = 1;
                        mk{end}.reldiv = 0;
                        %found(end+1) = 1;
                        bpms{end+1} = ptli;
                        foundk(end+1) = 1;
                        globpm(end+1,1:l) = NaN(1,l);
                        globpm(end,l) = ptli;
                        for i3 = 1:size(globpm,1)-1
                            if globpm(i3,l) == 0
                                globpm(i3,l) = globpm(i3,l-1);
                            end
                        end
                        %coord = [length(bpms),1];
                    end
                end
                
                active = zeros(1,length(mk));
                for i = 1:length(mk)
                    if ~mk{i}(1).active
                        continue
                    end
                    for i2 = 1:length(mk{i})
                        if isempty(mk{i}(i2).function) || ...
                                mk{i}(i2).timidx(end) < l 
                            continue
                        end
                        if isempty(mk{i}(i2).ref)
                            continue
                        end
                        ref = find([mk{i}.lvl] == mk{i}(i2).ref);
                        if isempty(mk{i}(ref).function) || ...
                                (mk{i}(ref).timidx(end) == l && ...
                                 mk{i}(ref).score(end) > mk{i}(i2).score(end))
                            continue
                        end
                        
                        reldiv = mk{i}(i2).reldiv;
                        same = find(mk{i}(ref).function(1,:) * reldiv < 0);
                        for i3 = 1:length(same)
                            refdiv = mk{i}(ref).function(1,same(i3));
                            otherlvl = mk{i}(ref).function(2,same(i3));
                            other = find([mk{i}.lvl] == otherlvl);
                            if isempty(other) || ...
                                    isempty(mk{i}(other).function) || ...
                                    abs(reldiv) >= abs(refdiv)
                                continue
                            end
                            intradiv = abs(refdiv/reldiv);
                            otherfunction = ...
                                find(mk{i}(other).function(1,:)...
                                     == -refdiv);
                            if isempty(otherfunction)
                                continue
                            end
                            if round(intradiv) == intradiv
                                mk{i}(ref)...
                                    .function(:,same(i3)) = ...
                                        [-reldiv; mk{i}(i2).lvl];
                                mk{i}(other).function(:,otherfunction)...
                                    = [intradiv; mk{i}(i2).lvl];
                                mk{i}(i2).function = ...
                                    [reldiv,-intradiv; ...
                                     mk{i}(ref).lvl, mk{i}(other).lvl];
                                break
                            elseif mk{i}(other).timidx(end) ~= l || ...
                                    mk{i}(other).score(end)...
                                        < mk{i}(i2).score(end)
                                mk{i}(i2).function = ...
                                    [reldiv; mk{i}(ref).lvl];
                                mk{i}(ref).function(:,end+1) = ...
                                    [-reldiv; mk{i}(i2).lvl];
                                mk{i}(other).function(:,otherfunction) = [];
                                break
                            end
                        end
                    end
                    
                    if l == 1 || isnan(globpm(i,l-1))
                        glo = 0;
                        sco = 0;
                        for i2 = 1:length(mk{i})
                            if mk{i}(i2).timidx(end) == l && ...
                                    ~isempty(mk{i}(i2).function)
                                sco2 = mk{i}(i2).score(end);
                                ind = mk{i}(i2).bpms(end) * mk{i}(i2).lvl;
                                glo = glo + ind * sco2;
                                sco = sco + sco2;
                            end
                        end
                        globpm(i,l) = glo / sco;
                        active(i) = 1;
                    else
                        glog = log2(globpm(i,l-1));
                        glodif = 0;
                        sco = 0;
                        mindif = Inf;
                        for i2 = 1:length(mk{i})
                            if mk{i}(i2).timidx(end) == l && ...
                                    ~isempty(mk{i}(i2).function) && ...
                                    mk{i}(i2).bpms(end) > 30
                                globpm2 = mk{i}(i2).bpms(end) ...
                                          * mk{i}(i2).lvl;
                                dif = glog - log2(globpm2);
                                sco2 = mk{i}(i2).score(end) / abs(dif);
                                glodif = glodif + dif * sco2;
                                sco = sco + sco2;
                                if abs(dif) < abs(mindif)
                                    mindif = dif;
                                end
                            end
                        end
                        if glodif
                            glodif = glodif / sco;
                            active(i) = 1;
                        else
                            active(i) = 0;
                        end
                        if 0 %abs(glodif) > abs(mindif)
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
                if ~isempty(find(active,1))
                    inactive = find(~active);
                    for i = 1:length(inactive)
                        mk{inactive(i)}(1).active = 0;
                    end
                end
                
                i = 0;
                while i < length(mk)
                    i = i + 1;
                    if mk{i}(1).locked
                        continue
                    end
                    if ~mk{i}(1).active
                        continue
                    end
                    
                    for i3 = 1:length(mk{i})
                        score1 = 0;
                        if mk{i}(i3).element
                            nbpms1 = globpm(i,l)/ mk{i}(i3).lvl;
                            lvl1 = mk{i}(i3).lvl;
                        end
                        if mk{i}(i3).score(end) > score1
                            score1 = mk{i}(i3).score(end);
                        end
                    end
                    
                    i3 = 1;
                    while i3 < i
                        if ~mk{i3}(1).active
                            i3 = i3 + 1;
                            continue
                        end
                        
                        included = 0;
                        score2 = 0;
                        for i2 = 1:length(mk{i3})
                            if mk{i3}(i2).element
                                nbpms2 = globpm(i3,l)/ mk{i3}(i2).lvl;
                                lvl2 = mk{i3}(i2).lvl;
                                if abs(60/nbpms1 - 60./nbpms2) < .01
                                    included = i2;
                                end
                            end
                            if mk{i3}(i2).score(end) > score2
                                score2 = mk{i3}(i2).score(end);
                            end
                        end
                        
                        if included
                            if score1 > score2
                                majo = i;
                                mino = i3;
                                ratio = lvl1 / lvl2;
                            else
                                majo = i3;
                                mino = i;
                                ratio = lvl2 / lvl1;
                            end
                            
                            for i2 = 1:length(mk{mino})
                                lvl = mk{mino}(i2).lvl * ratio;
                                i5 = find(lvl == [mk{majo}.lvl],1);
                                if isempty(i5)
                                    mk{majo}(end+1).lvl = lvl;
                                    mk{majo}(end).lastbpm = mk{mino}(i2).lastbpm;
                                    mk{majo}(end).bpms = mk{mino}(i2).bpms;
                                    mk{majo}(end).globpms = globpm(i,mk{mino}(i2).timidx) ...
                                        / mk{mino}(i2).lvl;
                                    mk{majo}(end).timidx = mk{mino}(i2).timidx;
                                    mk{majo}(end).score = mk{mino}(i2).score;
                                    mk{majo}(end).ref = mk{mino}(i2).ref * ratio;
                                    mk{majo}(end).reldiv = mk{mino}(i2).reldiv;
                                    bpms{majo}(end+1) = mk{mino}(i2).lastbpm;

                                end
                            end
                            
                            mk(mino) = [];
                            globpm(mino,:) = [];
                            bpms(mino) = [];
                            if score1 > score2
                                i3 = i3 - 1;
                            else
                                i = i - 1;
                                break
                            end
                        end
                        
                        i3 = i3 + 1;
                    end
                    
                    for i3 = 1:0 length(mk)
                        if i == i3 || ~mk{i3}(1).active
                            continue
                        end
                        included = 0;
                        for i2 = 1:length(mk{i})
                            if isempty(mk{i}(i2).function)
                                continue
                            end
                            
                            nbpms1 = globpm(i,l)/ mk{i}(i2).lvl;
                            nbpms2 = repmat(globpm(i3,l),[1,size(mk{i3},2)])...
                                                ./ [mk{i3}.lvl];
                            for i4 = 1:length(mk{i3})
                                if mk{i3}(i4).timidx(end) < l %|| ...
                                        %isempty(mk{i3}(i4).function)
                                    nbpms2(i4) = NaN;
                                end
                            end
                            dist = abs(60/nbpms1 - 60./nbpms2);
                            i4 = find(dist<.01,1);
                            if ~isempty(i4)
                                if isempty(mk{i3}(i4).function)
                                    continue
                                    mk{i}(1).locked = i3;
                                end
                                
                                included = 1;
                                ratio = mk{i3}(i4).lvl / mk{i}(i2).lvl;
                                for i4 = 1:length(mk{i})
                                    lvl = mk{i}(i4).lvl * ratio;
                                    i5 = find(lvl == [mk{i3}.lvl],1);
                                    if isempty(i5)
                                        mk{i3}(end+1).lvl = lvl;
                                        mk{i3}(end).lastbpm = mk{i}(i4).lastbpm;
                                        mk{i3}(end).bpms = mk{i}(i4).bpms;
                                        mk{i3}(end).globpms = globpm(i,mk{i}(i4).timidx) ...
                                            / mk{i}(i4).lvl;
                                        mk{i3}(end).timidx = mk{i}(i4).timidx;
                                        mk{i3}(end).score = mk{i}(i4).score;
                                        mk{i3}(end).ref = mk{i}(i4).ref * ratio;
                                        mk{i3}(end).reldiv = mk{i}(i4).reldiv;
                                        bpms{i3}(end+1) = mk{i}(i4).lastbpm;
                                    else
                                        new.timidx = [];
                                        new.bpms = [];
                                        new.globpms = [];
                                        new.score = [];
                                        for it = 1:l
                                            t2 = find(mk{i3}(i5).timidx...
                                                == it);
                                            t1 = find(mk{i}(i4).timidx...
                                                == it);
                                            if isempty(t1)
                                                if ~isempty(t2)
                                                    new.timidx(end+1) = it;
                                                    new.bpms(end+1) = ...
                                                        mk{i3}(i5).bpms(t2);
                                                    new.globpms(end+1) = ...
                                                        mk{i3}(i5).globpms(t2);
                                                    new.score(end+1) = ...
                                                        mk{i3}(i5).score(t2);
                                                end
                                            else
                                                if isempty(t2) || ...
                                                        mk{i3}(i5).score(t2) < ...
                                                            mk{i}(i4).score(t1)
                                                    test = 1;
                                                    for i6 = 1:length(mk{i3})
                                                        t3 = find(mk{i3}(i6).timidx == it);
                                                        if ~isempty(t3) && ...
                                                                abs(60./mk{i3}(i6).bpms(t3) ...
                                                                    - 60./mk{i}(i4).bpms(t1)) ...
                                                                    < .01
                                                            test = 0;
                                                            break
                                                        end
                                                    end
                                                else
                                                    test = 0;
                                                end
                                                if test
                                                    new.timidx(end+1) = it;
                                                    new.bpms(end+1) = ...
                                                        mk{i}(i4).bpms(t1);
                                                    new.globpms(end+1) = ...
                                                        mk{i}(i4).globpms(t1);
                                                    new.score(end+1) = ...
                                                        mk{i}(i4).score(t1);
                                                elseif ~isempty(t2)
                                                    new.timidx(end+1) = it;
                                                    new.bpms(end+1) = ...
                                                        mk{i3}(i5).bpms(t2);
                                                    new.globpms(end+1) = ...
                                                        mk{i3}(i5).globpms(t2);
                                                    new.score(end+1) = ...
                                                        mk{i3}(i5).score(t2);
                                                end
                                            end
                                        end
                                        mk{i3}(i5).timidx = new.timidx;
                                        mk{i3}(i5).bpms = new.bpms;
                                        mk{i3}(i5).globpms = new.globpms;
                                        mk{i3}(i5).score = new.score;
                                    end
                                end
                            end
                        end
                        if included
                            % meters{i} is completely included into
                            % meters{i2}
                            mk(i) = [];
                            globpm(i,:) = [];
                            i =  i - 1;
                            break
                        end
                    end
                end
            end
        end

        meters{j}{k} = mk;
        globpms{j}{k} = globpm;
    end
end

dm = purgedata(p);
m.autocor = dm;
m.globpm = globpms;
dm = set(dm,'Data',meters);
m = class(m,'mirmetre',mirdata(dm));


function bpm = getbpm(p,ptl)
if isa(p,'mirautocor') && not(get(p,'FreqDomain'))
    bpm = 60./ptl;
else
    bpm = ptl*60;
end