function varargout = mirentropy(x,varargin)
%   h = mirentropy(a) calculates the relative entropy of a.
%   (Cf. User's Manual.)
%   mirentropy(..., ?Center?) centers the input data before
%       transforming it into a probability distribution.

        center.key = 'Center';
        center.type = 'Boolean';
        center.default = 0;
    option.center = center;
    
        minrms.key = 'MinRMS';
        minrms.when = 'After';
        minrms.type = 'Numerical';
        minrms.default = .005;
    option.minrms = minrms;

specif.option = option;

varargout = mirfunction(@mirentropy,x,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if isamir(x,'miraudio')
    x = mirspectrum(x);
end
type = 'mirscalar';


function h = main(x,option,postoption)
if iscell(x)
    x = x{1};
end
m = get(x,'Data');
v = cell(1,length(m));
for h = 1:length(m)
    mh = m{h};
    if ~iscell(mh)
        mh = {mh};
    end
    v{h} = cell(1,length(mh));
    for k = 1:length(mh)
        mk = mh{k};
        mn = mk;
        if isa(x,'mirhisto') || isa(x,'mirscalar')
            mn = mn';
        end
        
        if option.center
            mn = center(mn);
        end
        
        % Negative data is trimmed:
        mn(mn<0) = 0;
        
        % Data is normalized such that the sum is equal to 1.
        mn = mn./repmat(sum(mn)+repmat(1e-12,...
                            [1 size(mn,2) size(mn,3) size(mn,4)]),...
                       [size(mn,1) 1 1 1]);
                   
        % Actual computation of entropy
        v{h}{k} = -sum(mn.*log(mn + 1e-12))./log(size(mn,1));
        
        if isa(x,'mirhisto') || isa(x,'mirscalar')
            v{h}{k} = v{h}{k}';
        end
    end
end
t = ['Entropy of ',get(x,'Title')];
h = mirscalar(x,'Data',v,'Title',t);
if (isa(x,'miraudio') || ...
    (isa(x,'mirspectrum')  && strcmpi(get(x,'Title'),'Spectrum')) || ...
    isa(x,'mircepstrum')) && ...
        isstruct(postoption) && ...
        isfield(postoption,'minrms') && postoption.minrms
    h = after(x,h,postoption.minrms);
end


function h = after(x,h,minrms)
r = mirrms(x,'Warning',0);
v = mircompute(@trim,get(h,'Data'),get(r,'Data'),minrms);
h = set(h,'Data',v);

    
function d = trim(d,r,minrms)
r = r/max(max(r));
pos = find(r<minrms);
for i = 1:length(pos)
    d(pos(i)) = NaN;
end
d = {d};