function varargout = mirmetroid(orig,varargin)
%   mc = mirmetroid(x) provides an estimation of dynamic metrical centroid.
%   [mc,ms] = mirmetroid(x) also provides an estimation of dynamic metrical
%   strength.
%
% When mirmetroid is used for academic research, please cite the following 
%   publication:
%   Lartillot, O., Cereghetti, D., Eliard, K., Trost, W. J., Rappaz, M.-A.,
%       Grandjean, D., "Estimating tempo and metrical features by tracking 
%       the whole metrical hierarchy", 3rd International Conference on 
%       Music & Emotion, Jyväskylä, 2013.
%
%   Optional arguments:
%       mirmetroid(..., ?Gate?) uses a simpler method, where the weights
%           are simply identified to the corresponding autocorrelation
%           scores, leading to possible abrupt changes in the metrical
%           centroid curve.

        gate.key = 'Gate';
        gate.type = 'Boolean';
        gate.default = 0;
    option.gate = gate;

        combine.key = 'Combine';
        combine.type = 'Boolean';
        combine.default = 1;
        combine.when = 'After';
    option.combine = combine;

specif.option = option;

varargout = mirfunction(@mirmetroid,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if ~isamir(x,'mirscalar') && ~isamir(x,'mirmetre')
    x = mirmetre(x);
end
type = {'mirscalar', 'mirscalar','mirmetre'};


function mo = main(m,option,postoption)
if iscell(m)
    m = m{1};
end
if isa(m,'mirscalar')
    mo = m;
else
    d = get(m,'Data');
    fp = get(m,'FramePos');
    [mo ms] = mircompute(@algo,d,fp,option.gate);
    mo = mirscalar(m,'Data',mo,'Title','Metrical Centroid','Unit','BPM');
    ms = mirscalar(m,'Data',ms,'Title','Metrical Strength');
end
[mo ms] = modif(mo,ms,postoption);
mo = {mo,ms,m};


function [m ms] = algo(d,fp,gate)
m = NaN(length(d),size(fp,2));
ms = zeros(length(d),size(fp,2));
for i = 1:length(d)
    timidx = [];
    for k = 1:length(d{i})
        timidx = union(timidx,d{i}(k).timidx);
    end
    lvl = [d{i}.lvl];
    [lvl ord] = sort(lvl);
    for j = 1:length(timidx)
        m(i,timidx(j)) = 0;
        ms(i,timidx(j)) = 0;
        t = NaN(1,length(lvl));
        for k = 1:length(lvl)
            tk = find(d{i}(ord(k)).timidx == timidx(j),1);
            if ~isempty(tk)
                t(k) = tk;
                if k == 1
                    submax = 0;
                else
                    sub = zeros(1,k-1);
                    mult = [];
                    for h = 1:k-1
                        if isempty(t(h)) || ...
                                isempty(find(d{i}(ord(h)).timidx == timidx(j),1))
                            continue
                        end
                        if ~mod(lvl(k),lvl(h))
                            sub(h) = d{i}(ord(h)).score(t(h));
                            mult(end+1) = h;
                        else
                            ismult = find(~mod(lvl(h),lvl(mult)),1);
                            if ~isempty(ismult)
                                h0 = mult(ismult);
                                div1 = round(lvl(h)/lvl(h0));
                                div2 = round(lvl(k)/lvl(h0));
                                if div2 - div1 == 1
                                    sub(h0) = min(sub(h0),...
                                                d{i}(ord(h)).score(t(h)));
                                end
                            end
                        end
                    end
                    submax = max(sub);
                end
                if gate
                    if d{i}(ord(k)).score(t(k)) > submax
                        scork = d{i}(ord(k)).score(t(k));
                    else
                        scork = 0;
                    end
                else
                    scork = max(0, d{i}(ord(k)).score(t(k)) - submax);
                end
                
                m(i,timidx(j)) = m(i,timidx(j)) ...
                                 + d{i}(ord(k)).bpms(t(k)) * scork;
                ms(i,timidx(j)) = ms(i,timidx(j)) + scork;
            end
        end
        m(i,timidx(j)) = m(i,timidx(j)) / ms(i,timidx(j));
    end
end


function [m ms] = modif(m,ms,option)
if option.combine
    dm = get(m,'Data');
    dms = get(ms,'Data');
    [dm dms] = mircompute(@combine,dm,dms);
    m = set(m,'Data',dm);
    ms = set(ms,'Data',dms);
end


function [m sms] = combine(m,ms)
m(isnan(m)) = 0;
sms = sum(ms,1);
m = sum(m.*ms,1)./sms;