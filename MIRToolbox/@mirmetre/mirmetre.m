function varargout = mirmetre(orig,varargin)
%   m = mirmetre(x) provides a detailed description of the hierarchical 
%       metrical structure by detecting periodicities from the onset 
%       detection curve and tracking a broad set of metrical levels.
%
% When mirmetre is used for academic research, please cite the following 
%   publication:
%   Lartillot, O., Cereghetti, D., Eliard, K., Trost, W. J., Rappaz, M.-A.,
%       Grandjean, D., "Estimating tempo and metrical features by tracking 
%       the whole metrical hierarchy", 3rd International Conference on 
%       Music & Emotion, Jyv?skyl?, 2013.
%
%   Optional arguments:
%       mirmetre(...,'Frame',l,h) orders a frame decomposition of window
%           length l (in seconds) and hop factor h, expressed relatively to
%           the window length. For instance h = 1 indicates no overlap.
%           Default values: l = 5 seconds and h = .05
%       mirmetre(..., ?Min?, mi) indicates the lowest periodicity taken 
%           into consideration, expressed in bpm. Default value: 24 bpm.
%       mirmetre(..., ?Max?, ma) specifies the highest periodicity taken 
%           into consideration, expressed in bpm. Default value: Inf, 
%           meaning that no limit is set a priori on the highest 
%           periodicity.
%       mirmetre(..., ?Contrast?, c) specifies the contrast factor for the 
%           peak picking. Default value: c = 0.05.
%       mirmetre(..., ?Threshold?, c) specifies the contrast factor for the
%           peak picking. Default value: c = 0.
%
%   [m,p] = mirmetre(...) also displays the autocorrelation function
%       leading to the tempo estimation, and shows in particular the
%       peaks corresponding to the tempo values.

        
        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        %frame.default = [0 0];
        frame.default = [5 .05];
    option.frame = frame;
    
        minres.key = 'MinRes';
        minres.type = 'Integer';
        minres.default = 10; %.1;
    option.minres = minres;
        
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
        ma.default = 1000; %500
    option.ma = ma;
    
        tol.key = 'Tolerance';
        tol.type = 'Integer';
        tol.default = .2;
    option.tol = tol;
    
        goto.key = {'Goto'};
        goto.type = 'Boolean';
        goto.default = 0;
    option.goto = goto;
        
specif.option = option;

varargout = mirfunction(@mirmetre,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if iscell(x)
    x = x{1};
end

if isamir(x,'mirmetre')
    return
end

if ~isamir(x,'mirautocor')
    if isamir(x,'mirenvelope')
        x = mironsets(x,'Frame',option.frame.length.val,...
                                option.frame.length.unit,...
                                option.frame.hop.val,...
                                option.frame.hop.unit);
    else
        if option.goto
            sg = 'Goto';
        else
            sg = 'Lartillot';
        end
        x = mironsets(x,'SmoothGate',sg,'MinRes',option.minres,...
                        'Detect',0,...
                        'Frame',option.frame.length.val,...
                                option.frame.length.unit,...
                                option.frame.hop.val,...
                                option.frame.hop.unit);
    end
end

x = mirautocor(x,'Min',60/option.ma,'Max',60/option.mi * 2,...
                 'NormalWindow',0);

x = mirpeaks(x,'Total',Inf,...
               'Threshold',option.thr,'Contrast',option.cthr,...
               'NoBegin','NoEnd',...
               'Normalize','Local','Order','Amplitude');
type = 'mirmetre';
    

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
            activestruct = []; %% Obsolete, activestruct can be entirely removed.
            globpm = [];
            for l = 1:size(ptk,2)       % For each successive frame
                %if ~mod(l,100)
                %    l
                %end
                
                ptl = getbpm(p,ptk{1,l,h}); % Peaks

                bpms = cell(1,length(mk));
                for i2 = 1:length(mk)
                    bpms{i2} = [mk{i2}.lastbpm];
                end
                
                foundk = zeros(1,length(mk));
                active = zeros(1,length(mk));
                new = 0; %zeros(1,length(mk));
                foundomin = zeros(1,length(mk));
                
                ampli = d{j}{k}(:,l,h);
                pos = ppp{j}{k}{l};
                
                %%                
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
                    
                    %thri = (1-(pv{j}{k}{l}(i) - mipv)/(mapv - mipv))^2/10 ...
                    %       + .1;
                    
                    score = ampli(pos(i));
                    found = zeros(1,length(mk));              % Is peak in metrical hierarchies?
                    
                    dist = inf(1,length(mk));
                    indx = nan(1,length(mk));
                    
                    i2 = 1;
                    while i2 <= length(mk)   % For each metrical hierarchy
                        if ~activestruct(i2) || ...
                                ~mk{i2}(1).active || isempty(bpms{i2})
                            i2 = i2+1;
                            continue
                        end
                        
                        globpmi2 = globpm(i2,end);
                        bpm2 = repmat(globpmi2, [1 length(mk{i2})])...
                               ./ [mk{i2}.lvl];
                        dist1 = abs(60/ptli - 60./bpm2);
                                                
                        dist2 = dist1;

                        if 0 % This short-cut has been toggled off, why?
                            for i3 = 1:length(mk{i2})
                                if isempty(mk{i2}(i3).function)
                                    dist1(i3) = NaN;
                                end
                            end

                            [disti1 indx1] = min(dist1);
                            if disti1 < thri && ...
                                    abs(log2(ptli / bpm2(indx1))) < .2
                                dist(i2) = disti1;
                                indx(i2) = indx1;
                                i2 = i2+1;
                                new = 0;
                                continue
                            end
                        end
                        
                        [disti2 indx2] = min(dist2);
                        
                        dist3 = NaN(1,length(mk{i2}));
                        for i3 = 1:length(mk{i2})
                            if mk{i2}(i3).timidx(end) == l
                                dist3(i3) = NaN; %unnecessary
                                continue
                            end
                            
                            t3 = find(mk{i2}(i3).timidx == l-1); %simplify
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
                        end
                        
                        if ~foundk(i2)
                            if abs(log2(ptli / bpm2(indx(i2)))) > .2 %.1
                                dist(i2) = Inf;
                                indx(i2) = 0;
                                new = 1; % New hierarchy (?)
                                i2 = i2+1;
                                continue
                            else
                                new = 0;
                            end
                        end
                        
                        i2 = i2+1;
                    end
                    
                    [unused order] = sort(dist);
                    for i2 = order
                        if ~activestruct(i2)
                            continue
                        end
                        
                        %if foundk(i2)
                        %    thri2 = thri;
                        %else
                        %    thri2 = min(thri,.1);
                        %end
                        thri2 = .07; %.01; %.07;
                        
                        if isnan(dist(i2)) || dist(i2) > thri2
                            continue
                        end
                        
                        % Continuing an existing metrical level.
                        
                        if mk{i2}(indx(i2)).timidx(end) ~= l
                            % Metrical level not extended yet.
                            
                            mk{i2}(indx(i2)).timidx(end+1) = l;
                            mk{i2}(indx(i2)).bpms(end+1) = ptl(i);
                            mk{i2}(indx(i2)).lastbpm = ptli;
                            mk{i2}(indx(i2)).score(end+1) = ...
                                d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h);

                            if foundk(i2)
                                active(i2) = 1;
                            else
                                % Metrical hierarchy not extended yet.
                                
                                if isempty(mk{i2}(indx(i2)).function)
                                    if isempty(find(foundk,1)) && ...
                                            mk{i2}(indx(i2)).score(end) > .15 %.3 %.15
                                        i3 = find([mk{i2}.lvl] == ...
                                                  mk{i2}(indx(i2)).ref,1);
                                        if ~isempty(mk{i2}(i3).function)
                                            mk{i2}(indx(i2)).function = ...
                                                [mk{i2}(indx(i2)).reldiv; ...
                                                 mk{i2}(indx(i2)).ref];
                                        end
                                    end
                                else
                                    active(i2) = 1;
                                end
                                foundk(i2) = 1;
                            end
                            
                            if ~foundomin(i2) && ...
                                    ~isempty(mk{i2}(indx(i2)).function)
                                % Global BPM determined using only the most
                                % dominant level.
                                foundomin(i2) = 1;
                                globpm(i2,l) = ptli * mk{i2}(indx(i2)).lvl;
                                for i3 = 1:size(globpm,1)
                                    if globpm(i3,l) == 0
                                        globpm(i3,l) = globpm(i3,l-1);
                                    end
                                end
                            end
                        end
                        found(i2) = foundk(i2);
                    end
                    
                    i2 = 1;
                    while i2 <= length(mk)
                        if ~activestruct(i2)
                            globpm(i2,l) = NaN;
                            for i3 = 1:size(globpm,1)
                                if globpm(i3,l) == 0
                                    globpm(i3,l) = globpm(i3,l-1);
                                end
                            end
                            i2 = i2+1;
                            continue
                        end
                        
                        %if found(i2)%%%%%%%%%%%%%%%%%%%%%%%%
                        %    i2 = i2+1;
                        %    continue
                        %end
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
                                    ...mk{i2}(ord(i3)).score(end) < .1 || ... %%%%% To toggle off sometimes? (cf. level 1/6 in the paper)
                                    ~isempty(mk{i2}(ord(i3)).complex)
                                i3 = i3-1;
                                continue
                            end
                            
                            if mk{i2}(ord(i3)).reldiv > 0 && ...
                                    isempty(mk{i2}(ord(i3)).function)
                                i3 = i3-1;
                                continue
                            end
                                                            
                            bpm3 = globpm(i2,end) / mk{i2}(ord(i3)).lvl;
                            
                            if 0 %abs(60/ptli - 60/bpm3) > dist(i2)  %%%%% To toggle off sometimes? (cf. level 1/6 in the paper)
                                i3 = i3-1;
                                continue
                            end
                                
                            rdiv = round(ptli / bpm3);
                            if rdiv == 1
                                i3 = i3-1;
                                continue
                            end
                            
                            lvl = mk{i2}(ord(i3)).lvl / rdiv;
                            if ~isempty(find(lvl == [mk{i2}.lvl],1))
                                i3 = i3-1;
                                continue
                            end
                                
                            if rdiv == 0 || rdiv > 8 || ...
                                    abs(60/ptli - ...
                                        60/(globpm(i2,end) / lvl))...
                                       > dist(i2)
                                i3 = i3-1;
                                continue
                            end
                                                        
                            %if 1 %~foundk(i2)
                                div = ptli ./ bpm3;
                            %else
                            %    div = [ptli2; ptli1] ./ bpm3;
                            %end
                            %if floor(div(1)) ~= floor(div(2))
                            %    newerr = 0;
                            %else
                                newerr = min(mod(div,1),1-mod(div,1));
                            %end
                            if ~foundk(i2)
                                thr = .02;
                            else
                                thr = option.tol;
                            end
                            
                            if newerr > thr
                                i3 = i3-1;
                                continue
                            end
                            
                            %if 0 && ~isempty(find([mk{i2}.lvl] > lvl & ...
                            %                 [mk{i2}.lvl] < mk{i2}(ord(i3)).lvl))
                            %    i3 = i3-1;
                            %    continue
                            %end
                            
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
                        if ~found(i2)
                            i3 = fo;
                            err = Inf;
                            while i3 <= length(orbpms)
                                if l - mk{i2}(ord(i3)).timidx(end) > 10 || ...
                                        ~isempty(mk{i2}(ord(i3)).complex)
                                    i3 = i3+1;
                                    continue
                                end

                                %if ...~foundk(i2) &&
                                %        isempty(mk{i2}(ord(i3)).function)
                                %    i3 = i3+1;
                                %    continue
                                %end

                                bpm3 = globpm(i2,end) / mk{i2}(ord(i3)).lvl;
                                                                
                                %if ~foundk(i2)
                                    div = bpm3 ./ [ptli ptli];
                                %else
                                %    div = bpm3 ./ [ptli2;ptli1];
                                %end
                                rdiv = round(bpm3 / ptli);
                                if rdiv == 1
                                    i3 = i3+1;
                                    continue
                                end
                                lvl = mk{i2}(ord(i3)).lvl * rdiv;
                                if ~isempty(find(lvl == [mk{i2}.lvl],1))
                                    i3 = i3+1;
                                    continue
                                end
                                if rdiv <= 1 || ...
                                        abs(60/ptli - ...
                                            60/(globpm(i2,end) / lvl)) ...
                                           > dist(i2)
                                    i3 = i3+1;
                                    continue
                                end

                                if floor(div(1)) < floor(div(2))
                                    newerr = 0;
                                else
                                    newerr = min(min(mod(div,1)),...
                                             min(1-mod(div,1)));
                                end
                                if ~foundk(i2)
                                    thr = .1; %.01;
                                else
                                    thr = option.tol;
                                end
                                if newerr < thr
                                    % Candidate level can be
                                    % integrated in this metrical
                                    % hierarchy
                                    if newerr < err
                                        if isempty(faster)
                                            faster.ref = ord(i3);
                                            faster.lvl = lvl;
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
                                            %ptli = mean([ptli1,ptli2]);
                                            err = newerr;
                                            %break
                                        elseif mk{i2}(ord(i3)).lvl * rdiv ...
                                                ~= faster.lvl
                                            faster.ref = ord(i3);
                                            faster.lvl = lvl;
                                            faster.rdiv = rdiv;
                                            %faster = [];
                                            %break
                                        end
                                    end
                                end
                                i3 = i3 + 1;
                            end
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
                        
                        active(i2) = 1;
                        found(i2) = 1;
                        new = 0;
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
                                        if round(intradiv) ~= intradiv
                                            continue
                                        end
                                        if ~isempty(mk{i2}(other).function)
                                            otherfunction = ...
                                                find(mk{i2}(other).function(1,:)...
                                                     == refdiv);
                                            if ~isempty(otherfunction)
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
                            if abs(reldiv) == 5 || abs(reldiv) > 6
                                mk{i2}(end).complex = 1;
                            end
                            
                            if reldiv < 0 && ~isempty(mk{i2}(ref).element)
                                mk{i2}(ref).element = [];
                                mk{i2}(end).element = 1;
                            end
                                                        
                            coord = [i2 length(mk{i2})];
                            bpms{i2}(end+1) = ptli;
                        else
                            dist0 = abs(mean([mk{i2}(l0).lastbpm]) - ptl(i));
                            [unused md] = min(dist0);
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
                        if ~foundk(i2)
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
                                        
                    if i == 1 && (new || ...
                                ~isempty(find(found == -1)) || ...
                                isempty(find(found))) && ...
                                d{j}{k}(ppp{j}{k}{1,l,h}(i),l,h) > .1 %.15)
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
                        mk{end}.complex = [];
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
                        %found = abs(found);
                        foundk = abs(foundk);
                        active(end+1) = 1;
                        %new(end+1) = 0;
                        foundomin(end+1) = 1;
                        activestruct(end+1) = 1;
                    end
                end
                
                %activestruct(~foundomin) = 0;
                
                %%
                for i = 1:length(mk)
                    %if ~activestruct(i) || ~mk{i}(1).active
                    %    globpm(i,l) = NaN;
                    %    for i2 = 1:length(mk{i})
                    %        if mk{i}(i2).timidx(end) == l
                    %            mk{i}(i2).globpms(end+1) = globpm(i,l) ...
                    %                                        / mk{i}(i2).lvl;
                    %        end
                    %    end
                    %    continue
                    %end
                    
                    %% Should we reactivate that test?
                    %if ~mk{i}(1).active
                    %    continue
                    %end
                    %%
                    
                    for i2 = 1:length(mk{i})
                        if isempty(mk{i}(i2).function) || ...
                                mk{i}(i2).timidx(end) < l 
                            continue
                        end
                        if isempty(mk{i}(i2).ref)
                            continue
                        end
                        ref = find([mk{i}.lvl] == mk{i}(i2).ref);
                        if isempty(ref)
                            continue
                        end
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
                                     == -refdiv,1);
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
                    
                    if l == 1 || isnan(globpm(i,l-1)) || ~globpm(i,l-1)
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
                    else
                        glog = log2(globpm(i,l-1));
                        glodif = 0;
                        sco = 0;
                        mindif = Inf;
                        for i2 = 1:length(mk{i})
                            if mk{i}(i2).timidx(end) == l && ...
                                    ~isempty(mk{i}(i2).function) && ...
                                    mk{i}(i2).score(end) > 0
                                %mk{i}(i2).bpms(end) > 30
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
                        end
                        if abs(glodif) > abs(mindif)
                            if glodif * mindif < 0
                                glodif = 0;
                            else
                                glodif = mindif;
                            end
                        end
                        if abs(glodif) > .05
                            glodif = .05 * sign(glodif);
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
                
                %%
                %if 1 %~isempty(find(active,1))
                    inactive = find(~active);
                    for i = 1:length(inactive)
                        mk{inactive(i)}(1).active = 0;
                    end
                %end
                
                %%
                i = 1;
                while i < length(mk)
                    i = i + 1;
                    if mk{i}(1).locked
                        continue
                    end
                    if ~mk{i}(1).active || ~activestruct(i)
                        continue
                    end
                    
                    score1 = 0;
                    for i3 = 1:length(mk{i})
                        if mk{i}(i3).element
                            nbpms1a = globpm(i,l)/ mk{i}(i3).lvl;
                            nbpms1b = mk{i}(i3).lastbpm;
                            element1 = i3;
                            %lvl1 = mk{i}(i3).lvl;
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
                            if ~isempty(mk{i3}(i2).element) || ...
                                    ~isempty(mk{i3}(i2).function)
                                nbpms2a = globpm(i3,l)/ mk{i3}(i2).lvl;
                                nbpms2b = mk{i3}(i2).lastbpm;
                                %lvl2 = mk{i3}(i2).lvl;
                                if abs(60/nbpms1a - 60./nbpms2a) < .01 || ...
                                        0 %abs(60/nbpms1b - 60./nbpms2b) < .01
                                    included = i2;
                                    element2 = i2;
                                end
                            end
                            if mk{i3}(i2).score(end) > score2
                                score2 = mk{i3}(i2).score(end);
                            end
                        end
                        
                        if included
                            if score1 > score2
                                majo = [i element1];
                                mino = i3;
                                if length(mk{i3}(1).timidx) > 1
                                    mk{i3}(1).active = 0;
                                    i3 = i3 + 1;
                                    break;
                                end
                            else
                                majo = [i3 element2];
                                mino = i;
                                if length(mk{i}(1).timidx) > 1
                                    mk{i}(1).active = 0;
                                    i3 = i3 + 1;
                                    break;
                                end
                            end
                            
                            for i2 = 1:length(mk{mino})
                                div = round(mk{majo(1)}(majo(2)).lastbpm ...
                                            / mk{mino}(i2).lastbpm);
                                lvl =  div * mk{majo(1)}(majo(2)).lvl;
                                i5 = find(lvl == [mk{majo(1)}.lvl],1);
                                if isempty(i5)
                                    mk{majo(1)}(end+1).lvl = lvl;
                                    mk{majo(1)}(end).lastbpm = mk{mino}(i2).lastbpm;
                                    mk{majo(1)}(end).bpms = mk{mino}(i2).bpms;
                                    mk{majo(1)}(end).globpms = ...
                                        globpm(majo(1), mk{mino}(i2).timidx)...
                                        / lvl;
                                    mk{majo(1)}(end).timidx = mk{mino}(i2).timidx;
                                    mk{majo(1)}(end).score = mk{mino}(i2).score;
                                    mk{majo(1)}(end).ref = mk{majo(1)}(majo(2)).lvl;
                                    mk{majo(1)}(end).reldiv = div;
                                    bpms{majo(1)}(end+1) = mk{mino}(i2).lastbpm;

                                end
                            end
                            
                            if score1 > score2
                                mk(i) = mk(i3);
                                globpm(i) = globpm(i3);
                                bpms(i) = bpms(i3);
                            end
                            
                            mk(i) = [];
                            globpm(i,:) = [];
                            bpms(i) = [];
                            active(i) = [];
                            %new(i) = [];
                            
                            i = i - 1;
                            break
                        end
                        
                        i3 = i3 + 1;
                    end
                    
                    %for i3 = 1:0 length(mk)
                    %    if i == i3 || ~mk{i3}(1).active
                    %        continue
                    %    end
                    %    included = 0;
                    %    for i2 = 1:length(mk{i})
                    %        if isempty(mk{i}(i2).function)
                    %            continue
                    %        end
                    %        
                    %        nbpms1 = globpm(i,l)/ mk{i}(i2).lvl;
                    %        nbpms2 = repmat(globpm(i3,l),[1,size(mk{i3},2)])...
                    %                            ./ [mk{i3}.lvl];
                    %        for i4 = 1:length(mk{i3})
                    %            if mk{i3}(i4).timidx(end) < l %|| ...
                    %                    %isempty(mk{i3}(i4).function)
                    %                nbpms2(i4) = NaN;
                    %            end
                    %        end
                    %        dist = abs(60/nbpms1 - 60./nbpms2);
                    %        i4 = find(dist<.01,1);
                    %        if ~isempty(i4)
                    %            if isempty(mk{i3}(i4).function)
                    %                continue
                    %                mk{i}(1).locked = i3;
                    %            end
                    %            
                    %            included = 1;
                    %            ratio = mk{i3}(i4).lvl / mk{i}(i2).lvl;
                    %            for i4 = 1:length(mk{i})
                    %                lvl = mk{i}(i4).lvl * ratio;
                    %                i5 = find(lvl == [mk{i3}.lvl],1);
                    %                if isempty(i5)
                    %                    mk{i3}(end+1).lvl = lvl;
                    %                    mk{i3}(end).lastbpm = mk{i}(i4).lastbpm;
                    %                    mk{i3}(end).bpms = mk{i}(i4).bpms;
                    %                    mk{i3}(end).globpms = globpm(i,mk{i}(i4).timidx) ...
                    %                        / mk{i}(i4).lvl;
                    %                    mk{i3}(end).timidx = mk{i}(i4).timidx;
                    %                    mk{i3}(end).score = mk{i}(i4).score;
                    %                    mk{i3}(end).ref = mk{i}(i4).ref * ratio;
                    %                    mk{i3}(end).reldiv = mk{i}(i4).reldiv;
                    %                    bpms{i3}(end+1) = mk{i}(i4).lastbpm;
                    %                else
                    %                    new.timidx = [];
                    %                    new.bpms = [];
                    %                    new.globpms = [];
                    %                    new.score = [];
                    %                    for it = 1:l
                    %                        t2 = find(mk{i3}(i5).timidx...
                    %                            == it);
                    %                        t1 = find(mk{i}(i4).timidx...
                    %                            == it);
                    %                        if isempty(t1)
                    %                            if ~isempty(t2)
                    %                                new.timidx(end+1) = it;
                    %                                new.bpms(end+1) = ...
                    %                                    mk{i3}(i5).bpms(t2);
                    %                                new.globpms(end+1) = ...
                    %                                    mk{i3}(i5).globpms(t2);
                    %                                new.score(end+1) = ...
                    %                                    mk{i3}(i5).score(t2);
                    %                            end
                    %                        else
                    %                            if isempty(t2) || ...
                    %                                    mk{i3}(i5).score(t2) < ...
                    %                                        mk{i}(i4).score(t1)
                    %                                test = 1;
                    %                                for i6 = 1:length(mk{i3})
                    %                                    t3 = find(mk{i3}(i6).timidx == it);
                    %                                    if ~isempty(t3) && ...
                    %                                            abs(60./mk{i3}(i6).bpms(t3) ...
                    %                                                - 60./mk{i}(i4).bpms(t1)) ...
                    %                                                < .01
                    %                                        test = 0;
                    %                                        break
                    %                                    end
                    %                                end
                    %                            else
                    %                                test = 0;
                    %                            end
                    %                            if test
                    %                                new.timidx(end+1) = it;
                    %                                new.bpms(end+1) = ...
                    %                                    mk{i}(i4).bpms(t1);
                    %                                new.globpms(end+1) = ...
                    %                                    mk{i}(i4).globpms(t1);
                    %                                new.score(end+1) = ...
                    %                                    mk{i}(i4).score(t1);
                    %                            elseif ~isempty(t2)
                    %                                new.timidx(end+1) = it;
                    %                                new.bpms(end+1) = ...
                    %                                    mk{i3}(i5).bpms(t2);
                    %                                new.globpms(end+1) = ...
                    %                                    mk{i3}(i5).globpms(t2);
                    %                                new.score(end+1) = ...
                    %                                    mk{i3}(i5).score(t2);
                    %                            end
                    %                        end
                    %                    end
                    %                    mk{i3}(i5).timidx = new.timidx;
                    %                    mk{i3}(i5).bpms = new.bpms;
                    %                    mk{i3}(i5).globpms = new.globpms;
                    %                    mk{i3}(i5).score = new.score;
                    %                end
                    %            end
                    %        end
                    %    end
                    %    if included
                    %        % meters{i} is completely included into
                    %        % meters{i2}
                    %        mk(i) = [];
                    %        globpm(i,:) = [];
                    %        i =  i - 1;
                    %        break
                    %    end
                    %end
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