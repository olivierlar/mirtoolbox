function varargout = mirtonalcentroid(orig,varargin)
%   c = mirtonalcentroid(x) calculates the 6-dimensional tonal centroid
%       vector from the chromagram. 
%   It corresponds to a projection of the chords along circles of fifths, 
%       of minor thirds, and of major thirds.
%   [c ch] = mirtonalcentroid(x) also returns the intermediate chromagram.
%
% C. A. Harte and M. B. Sandler, Detecting harmonic change in musical
%   audio, in Proceedings of Audio and Music Computing for Multimedia
%   Workshop, Santa Barbara, CA, 2006. 

        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [.743 .1];
    option.frame = frame;

specif.option = option;

varargout = mirfunction(@mirtonalcentroid,orig,varargin,nargout,specif,@init,@main);


function [c type] = init(orig,option)
if option.frame.length.val
    c = mirchromagram(orig,'Frame',option.frame.length.val,...
                                   option.frame.length.unit,...
                                   option.frame.hop.val,...
                                   option.frame.hop.unit,...
                                   option.frame.phase.val,...
                                   option.frame.phase.unit,...
                                   option.frame.phase.atend);
else
    c = mirchromagram(orig);
end
type = 'mirtonalcentroid';


function tc = main(ch,option,postoption)
if iscell(ch)
    ch = ch{1};
end
if isa(ch,'mirtonalcentroid')
    tc = orig;
    ch = [];
else
    x1 = sin(pi*7*(0:11)/6)';
    y1 = cos(pi*7*(0:11)/6)';
    % minor thirds circle
    x2 = sin(pi*3*(0:11)/2)';
    y2 = cos(pi*3*(0:11)/2)';
    % major thirds circle
    x3 = 0.5 * sin(pi*2*(0:11)/3)';
    y3 = 0.5 * cos(pi*2*(0:11)/3)';
    c = [x1 y1 x2 y2 x3 y3];
    c = c';
    tc = class(struct,'mirtonalcentroid',mirdata(ch));
    tc = purgedata(tc);
    tc = set(tc,'Title','Tonal centroid','Abs','dimensions','Ord','position');
    m = get(ch,'Magnitude');
    %disp('Computing tonal centroid...')
    n = cell(1,length(m));  % The final structured list of magnitudes.
    d = cell(1,length(m));  % The final structured list of centroid dimensions.
    for i = 1:length(m)
        mi = m{i};
        if not(iscell(mi))
            mi = {mi};
        end
        ni = cell(1,length(mi));    % The list of magnitudes.
        di = cell(1,length(mi));    % The list of centroid dimensions.
        for j = 1:length(mi)
            mj = mi{j};
            ni{j} = zeros(6,size(mj,2),size(mi,3));
            for k = 1:size(mj,3)
                ni{j}(:,:,k) = c * mj(:,:,k);
            end
            di{j} = repmat((1:6)',[1,size(mj,2),size(mi,3)]);
        end
        n{i} = ni;
        d{i} = di;
    end
    tc = set(tc,'Positions',n,'Dimensions',d);
end
tc = {tc,ch};