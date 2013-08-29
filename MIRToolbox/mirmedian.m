function varargout = mirmedian(f,varargin)
% m = mirmedian(f) returns the median along frames of the feature f
%
%   f can be a structure array composed of features. In this case,
%       m will be structured the same way.

if isa(f,'mirstruct')
    data = get(f,'Data');
    for fi = 1:length(data)
        data{fi} = mirmedian(data{fi});
    end
    varargout = {set(f,'Data',data)};
elseif isstruct(f)
    fields = fieldnames(f);
    for i = 1:length(fields)
        field = fields{i};
        stat.(field) = mirmedian(f.(field));
    end
    varargout = {stat};
else
    specif.nochunk = 1;
    varargout = mirfunction(@mirmedian,f,varargin,nargout,specif,@init,@main);
end


function [x type] = init(x,option)
type = '';


function m = main(f,option,postoption)
if iscell(f)
    f = f{1};
end
if isa(f,'mirhisto')
    warning('WARNING IN MIRmedian: histograms are not taken into consideration yet.')
    m = struct;
    return
end
fp = get(f,'FramePos');
ti = get(f,'Title');
if 0 %get(f,'Peaks')
    if not(isempty(get(f,'PeakPrecisePos')))
        stat = addstat(struct,get(f,'PeakPrecisePos'),fp,'PeakPos');
        stat = addstat(stat,get(f,'PeakPreciseVal'),fp,'PeakMag');
    else
        stat = addstat(struct,get(f,'PeakPosUnit'),fp,'PeakPos');
        stat = addstat(stat,get(f,'PeakVal'),fp,'PeakMag');
    end
else
    d = get(f,'Data');
end
l = length(d);
m = cell(1,l);
for i = 1:l
    dd = d{i};
    if iscell(dd)
        m{i} = cell(1,length(dd));
        fpi = cell(1,length(dd));
        for j = 1:length(dd)
            ddj = dd{j};
            if iscell(ddj)
                ddj = uncell(ddj);
            end
            for k = 1:size(ddj,1)
                ddk = ddj(k,:);
                ddk(isnan(ddk)) = [];
                m{i}{j}(k,1) = median(ddk,2);
            end
            fpi{j} = [fp{i}{j}(1);fp{i}{j}(end)];
        end
        fp{i} = fpi;
    elseif size(dd,2) < 2
        nonan = find(not(isnan(dd)));
        dn = dd(nonan);
        m{i}{1} = median(dn,2);
    else
        %diffp = fp{i}{1}(1,2:end) - fp{i}{1}(1,1:end-1);
        %if round((diffp(2:end)-diffp(1:end-1))*1000)
            % Not regular sampling (in mirattacktime for instance)
        %    framesampling = NaN;
        %else
        %    framesampling = fp{i}{1}(1,2)-fp{i}{1}(1,1);
        %end
        dd = median(dd,4);
        m{i} = {NaN(size(dd,1),1,size(dd,3))};
        for k = 1:size(dd,1)
            for l = 1:size(dd,3)
                dk = dd(k,:,l);
                nonan = find(not(isnan(dk)));
                if not(isempty(nonan))
                    dn = dk(nonan);
                    m{i}{1}(k,1,l) = median(dn,2);
                end
            end
        end
    end
end
m = mirscalar(f,'Data',m,'Title',['Median of ',ti],'FramePos',fp);