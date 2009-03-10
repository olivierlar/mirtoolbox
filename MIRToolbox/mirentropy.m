function varargout = mirentropy(x,varargin)
%   h = mirentropy(a) calculates the relative entropy of a.
% 
%© Part of the MIDI Toolbox, Copyright © 2004, University of Jyvaskyla, Finland
% Reworked for MIR toolbox


varargout = mirfunction(@mirentropy,x,varargin,nargout,struct,@init,@main);


function [x type] = init(x,option)
if not(isamir(x,'mirdata'))
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
    v{h} = cell(1,length(m{h}));
    for k = 1:length(m{h})
        mk = m{h}{k};
        mn = mk;
        if isa(x,'mirhisto') || isa(x,'mirscalar')
            mn = mn';
        end
        mn = center(mn);
        mn(mn<0) = 0;
        %if min(min(min(mn)))<0
        %    mn = mn-repmat(min(mn),[size(mn,1) 1 1]);
        %end
        mn = mn./repmat(sum(mn)+repmat(1e-12,...
                            [1 size(mn,2) size(mn,3) size(mn,4)]),...
                       [size(mn,1) 1 1 1]);
        lgn = log(mn + 1e-12);
        v{h}{k} = -sum(mn.*lgn)./log(size(mn,1));
        if isa(x,'mirhisto') || isa(x,'mirscalar')
            v{h}{k} = v{h}{k}';
        end
    end
end
t = ['Entropy of ',get(x,'Title')];
h = mirscalar(x,'Data',v,'Title',t);