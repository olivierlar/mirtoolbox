function varargout = mirfluctuation(orig,varargin)
%   f = mirfluctuation(x) calculates the fluctuation strength, indicating the
%       rhythmic periodicities along the different channels.
%   Optional arguments:
%       mirfluctuation(...,'MinRes',mr) specifies the minimal frequency
%           resolution of the resulting spectral decomposition, in Hz.
%               Default: mr = .01 Hz
%       mirfluctuation(...,'Summary') returns the summary of the fluctuation,
%           i.e., the summation along the critical bands.
%
% E. Pampalk, A. Rauber, D. Merkl, "Content-based Organization and 
% Visualization of Music Archives", 

        sum.key = 'Summary';
        sum.type = 'Boolean';
        sum.default = 0;
    option.sum = sum;

        mr.key = 'MinRes';
        mr.type = 'Integer';
        mr.default = .01;
    option.mr = mr;
    
        inframe.key = 'InnerFrame';
        inframe.type = 'Integer';
        inframe.number = 2;
        inframe.default = [.023 .5];
    option.inframe = inframe;
    
        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [100 50];
    option.frame = frame;

specif.option = option;
     
varargout = mirfunction(@mirfluctuation,orig,varargin,nargout,specif,@init,@main);


function [s type] = init(x,option)
if iscell(x)
    x = x{1};
end
if isamir(x,'miraudio') && not(isframed(x))
    x = mirframe(x,option.inframe(1),option.inframe(2));
end
s = mirspectrum(x,'Power','Terhardt','Bark','dB','Mask');
type = 'mirspectrum';


function f = main(x,option,postoption)
d = get(x,'Data');
fp = get(x,'FramePos');
fl = option.frame.length.val;
fh = option.frame.hop.val;
if ~fl
    f = mirspectrum(x,'AlongBands','Max',10,'Window',0,'NormalLength',...
                      'Resonance','Fluctuation','MinRes',option.mr);
else
    for i = 1:length(d)
        for j = 1:length(d{i})
            if ~fl
                n = 1;
            else
                n = floor((size(d{i}{j},2)-fl)/fh)+1;  % Number of frames
            end
            for k = 1:n   % For each frame, ...
                if ~fl
                    st = 1;
                    stend = size(d{i}{j},2);
                else
                    st = (k-1)*fh+1;
                    stend = st+fl-1;
                end
                dk = d{i}{j}(:,st:stend,:);
                fpk = fp{i}{j}(:,st:stend);
                x2 = set(x,'Data',{{dk}},'FramePos',{{fpk}});
                fk = mirspectrum(x2,'AlongBands','Max',10,'Window',0,...
                                    'NormalLength',...
                                    'Resonance','Fluctuation',...
                                    'MinRes',option.mr);
                df{i}{j}(:,k,:) = mirgetdata(fk);
                fp2{i}{j}(:,k) = [fpk(1);fpk(end)];
            end
        end
    end
    f = set(fk,'Data',df,'FramePos',fp2);
end

if option.sum
    f = mirsummary(f);
end
f = set(f,'Title','Fluctuation');