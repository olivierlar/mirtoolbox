function varargout = mirstat(f,varargin)
% stat = mirstat(f) returns basic statistics of the feature f as a
%   structure array stat such that:
%       stat.Mean is the mean of the feature
%       stat.Std is the standard deviation
%       stat.Slope is the slope
%       stat.PeriodFreq is the frequency (in Hz) of the main periodicity
%       stat.PeriodAmp is the amplitude of the main periodicity
%       stat.PeriodEntropy is the entropy of the periodicity curve (which 
%           is identified to the autocorrelation function of f)
%
%   f can be itself a structure array composed of features. In this case,
%       stat will be structured the same way.
%
% mirstat does not work for multi-channels objects.

if isa(f,'mirstruct')
    data = get(f,'Data');
    for fi = 1:length(data)
        data{fi} = mirstat(data{fi});
    end
    varargout = {set(f,'Data',data)};
elseif isstruct(f)
    fields = fieldnames(f);
    for i = 1:length(fields)
        field = fields{i};
        stat.(field) = mirstat(f.(field));
    end
    varargout = {stat};
else
    specif.nochunk = 1;
    varargout = mirfunction(@mirstat,f,varargin,nargout,specif,@init,@main);
end


function [x type] = init(x,option)
type = '';


function stat = main(f,option,postoption)
if iscell(f)
    f = f{1};
end
if isa(f,'mirhisto')
    warning('WARNING IN MIRSTAT: histograms are not taken into consideration yet.')
    stat = struct;
    return
end
fp = get(f,'FramePos');
ti = get(f,'Title');
if haspeaks(f)
    ppp = get(f,'PeakPrecisePos');
    if not(isempty(ppp)) && not(isempty(ppp{1}))
        stat = addstat(struct,ppp,fp,'PeakPos');
        stat = addstat(stat,get(f,'PeakPreciseVal'),fp,'PeakMag');
    else
        if isa(f,'mirkeystrength') || isa(f,'mirchromagram')
            stat = addstat(struct,get(f,'PeakPos'),fp,'PeakPos');
        else
            stat = addstat(struct,get(f,'PeakPosUnit'),fp,'PeakPos');
        end
        stat = addstat(stat,get(f,'PeakVal'),fp,'PeakMag');
    end
else
    stat = addstat(struct,get(f,'Data'),fp,'');
end


function stat = addstat(stat,d,fp,field)
l = length(d);
anyframe = 0;
for i = 1:l
    if not(iscell(d{i}))
        d{i} = {d{i}};
    end
    for k = 1:length(d{i})
        dd = d{i}{k};
        if iscell(dd)
            if length(dd)>1
                anyframe = 1;
            end
            dn = zeros(1,length(dd));
            for j = 1:length(dd)
                if isempty(dd{j})
                    dn(j) = NaN;
                else
                    dn(j) = dd{j}(1);
                end
            end
            [m{i}{k},s{i}{k},sl{i}{k},pa{i}{k},pf{i}{k},pe{i}{k}] ...
                = computestat(dn,fp{i}{k});
        elseif size(dd,2) < 2
            nonan = find(not(isnan(dd)));
            dn = dd(nonan);
            if isempty(dn)
                m{i}{k} = NaN;
            else
                m{i}{k} = mean(dn,2);
            end
            s{i}{k} = NaN;
            sl{i}{k} = NaN;
            pa{i}{k} = NaN;
            pf{i}{k} = NaN;
            pe{i}{k} = NaN;
        else
            anyframe = 1;
            dd = mean(mean(dd,3),4);
            m{i}{k} = NaN(size(dd,1),1);
            s{i}{k} = NaN(size(dd,1),1);
            sl{i}{k} = NaN(size(dd,1),1);
            pa{i}{k} = NaN(size(dd,1),1);
            pf{i}{k} = NaN(size(dd,1),1);
            pe{i}{k} = NaN(size(dd,1),1);
            for j = 1:size(dd,1)
                [m{i}{k}(j),s{i}{k}(j),sl{i}{k}(j),...
                 pa{i}{k}(j),pf{i}{k}(j),pe{i}{k}(j)] = ...
                    computestat(dd(j,:),fp{i}{k});
            end
        end
    end
end
if anyframe
    fields = {'Mean','Std','Slope','PeriodFreq','PeriodAmp','PeriodEntropy'};
    stats = {m,s,sl,pf,pa,pe};   
else
    fields = {'Mean'};
    stats = {m};   
end
for i = 1:length(stats)
    data = stats{i};
    data = uncell(data,NaN);
    stat.(strcat(field,fields{i})) = data;
end


function [m,s,sl,pa,pf,pe] = computestat(d,fp)
m = NaN;
s = NaN;
sl = NaN;
pa = NaN;
pf = NaN;
pe = NaN;
diffp = fp(1,2:end) - fp(1,1:end-1);
if isempty(diffp) || sum(round((diffp(2:end)-diffp(1:end-1))*1000))
    % Not regular sampling (in mirattacktime for instance)
    framesampling = NaN;
else
    framesampling = fp(1,2)-fp(1,1);
end
nonan = find(not(isnan(d)));
if not(isempty(nonan))
    dn = d(nonan);
    m = mean(dn,2);
    s = std(dn,0,2);
    if not(isnan(s))
        if s
            dk = (dn-m)/s;
            tk = linspace(0,1,size(d,2));
            sl = dk(:)'/tk(nonan);
        elseif size(s,2) == 1
            s = NaN;
        end
    end
    if length(dn)>1
        cor = xcorr(dn',dn','coeff');
        cor = cor(ceil(length(cor)/2):end);
        % let's zero the first descending slope of the
        % autocorrelation function
        firstmin = find(cor(2:end)>cor(1:end-1));
        if not(isempty(firstmin) || isnan(framesampling))
            cor2 = cor;
            cor2(1:firstmin(1)) = 0;
            [pa,pfk] = max(cor2);
            if pfk > 1
                pf = 1/(pfk-1)/framesampling;
            end
        end
        cor = cor-min(cor);
        cor = cor/sum(cor);
        pe = -sum(cor.*log(cor+1e-12))./log(length(cor));
    end
end