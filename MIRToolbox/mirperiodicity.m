function [r,s,a,f] = mirperiodicity(x,varargin)
%   r = mirperiodicity(x) estimates the periodicity.
%           <<<< UNDER CONSTRUCTION!!!! >>>>
%   Optional arguments:
%       r = mirperiodicity(...,'Frame',l,h): orders a frame decomposition of
%           the audio input of window length l (in seconds) and hop factor
%           h, expressed relatively to the window length.
%           For instance, h = 1 indicates no overlap.
%           Default values: l = 3 seconds and h = .1
%       r = mirperiodicity(...,'Min',m): discards periodicities lower than
%           m Hz. Default value: 0.
%       r = mirperiodicity(...,'Max',m): discards periodicities higher than
%           m Hz. Default value: Inf.
%   [r,s] = mirperiodicity(x) also returns the inter-beat spectrum used for
%       the computation.
%   [r,s,a] = mirperiodicity(x) also returns the beat autocorrelation.

s = {};
a = {};
f = {};
para = scanargin(varargin);
if isa(x,'mirautocor') 
    a = x;
else
    [t,a,f] = mirtempo(x,'Frame',para.fl,para.fh);
end

if 0
    s = mirspectrum(a,'Min',para.mi,'Max',para.ma,'ZeroPad')
    m = get(s,'Magnitude');
    ph = get(s,'Phase');
    %p = peaks(s,'Total',1);
    %pp = peakpos(p);
    rc = cell(1,length(m));
    for l = 1:length(m)
        rc = cell(1,length(m{l}));
        for k = 1:length(m{l})
            %ppk = pp{k}{1};
            mk = m{l}{k};
            %phk = abs(pi - abs(ph{l}{k}));
            mpk = mk; %.*phk;
            [pv,pp] = max(mpk);
            rck = zeros(1,size(mk,2),size(mk,3));
            for j = 1:size(mk,3)
                for i = 1:size(mk,2)
                    range = max(1,pp(1,i,j)-1):min(size(mk,1),pp(1,i,j)+1);
                    rck(1,i,j) = sum(mk(range,i,j))/sum(mk(:,i,j));
                end
            end
            rc{l}{k} = rck;
        end    
    end
elseif 0
    m = get(a,'Data');
    rc = cell(1,length(m));
    for l = 1:length(m)
        rc = cell(1,length(m{l}));
        for k = 1:length(m{l})
            mk = m{l}{k};
            l2 = floor(size(mk,1)/2);    % Number of autocorrelation coefficients divided by scaling factor
            t2 = l2*2;                  % The new size of the autocorrelation vector 
                                % is equal to the size of the scaled autocorrelation  
            mk2 = interp1(mk(1:t2,:,:),(1/2:1/2:l2)');
            mk2(find(isnan(mk2))) = 0;    % All the NaN values are changed into 0
            rck = zeros(1,size(mk,2),size(mk,3));
            for j = 1:size(mk,3)
                for i = 1:size(mk,2)
                    rck(1,i,j) = mk(1:t2,i,j)'*mk2(:,i,j)/t2;
                    figure
                    plot(mk(1:t2,i,j))
                    hold on
                    plot(mk2(:,i,j),'r')
                end
            end
            rc{l}{k} = rck;
        end
    end
else
    %a = mirautocor(a,'Min',22,'Hz');
    m = get(a,'Data');
    rc = cell(1,length(m));
    for l = 1:length(m)
        rc = cell(1,length(m{l}));
        for k = 1:length(m{l})
            mk = m{l}{k};
            rck = zeros(1,size(mk,2),size(mk,3));
            for j = 1:size(mk,3)
                for i = 1:size(mk,2)
                    ind = find(mk(:,i,j)<0);
                    if isempty(ind) % no zero-crossings
                        rck(1,i,j) = 0;
                    else
                        rck(1,i,j) = max(mk(ind(1):end,i,j))/mk(1,i,j); % search maximum from values after the first zro-crossing of ac3
                    end
                end
            end
            rc{l}{k} = rck;
        end
    end
end
r = mirscalar(a,'Data',rc,'Title','Periodicity');


function para = scanargin(v)
para.fl = 0;
para.fh = 0;
para.mi = 0;
para.ma = Inf;
if exist('v')
    i = 1;
    while i <= length(v)
        arg = v{i};
        if strcmpi(arg,'Frame')
            fl = 3;
            fh = .1;
            if length(v)>i && isnumeric(v{i+1})
                i = i+1;
                fl = v{i};
            end
            if length(v)>i && isnumeric(v{i+1})
                i = i+1;
                fh = v{i};
            end
        elseif strcmpi(arg,'Min')
            i = i+1;
            mi = v{i};
        elseif strcmpi(arg,'Max')
            i = i+1;
            ma = v{i};
        else
            error('ERROR IN MIRPERIODICITY: Syntax error. See help mirperiodicity.');
        end    
        i = i+1;
    end
end