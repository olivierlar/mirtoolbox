function varargout = mirfluctuation(orig,varargin)
%   f = mirfluctuation(x) calculates the fluctuation strength, indicating
%       the rhythmic periodicities along the different channels.
%   Optional arguments:
%       mirfluctuation(...,'MinRes',mr) specifies the minimal frequency
%           resolution of the resulting spectral decomposition, in Hz.
%               Default: mr = .01 Hz
%       mirfluctuation(...,'Summary') returns the summary of the
%           fluctuation, i.e., the summation along the critical bands.
%       mirfluctuation(..., 'InnerFrame', l, r) specifies the spectrogram 
%           frame length l (in second), and, optionally, the frame rate r 
%           (in Hertz), with by default a frame length of 23 ms and a frame
%           rate of 80 Hz.
%       mirfluctuation(..., 'Frame', l, r) computes fluctuation using a 
%           window moving along the spectrogram, whose length l (in second)
%           and frame rate r (in Hertz) can be specified as well, with by
%           default a frame length of 1 s and a frame rate of 10 Hz.
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
    
        max.key = 'Max';
        max.type = 'Integer';
        max.default = 10;
    option.max = max;

        band.type = 'String';
        band.choice = {'Mel','Bark'};
        band.default = 'Bark';
    option.band = band;
    
        inframe.key = 'InnerFrame';
        inframe.type = 'Integer';
        inframe.number = 2;
        inframe.default = [.023 80];
    option.inframe = inframe;
    
        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [1 10];
    option.frame = frame;

specif.option = option;
specif.nochunk = 1;
     
varargout = mirfunction(@mirfluctuation,orig,varargin,nargout,specif,@init,@main);


function [s type] = init(x,option)
if iscell(x)
    x = x{1};
end
if option.inframe(2) < option.max * 2
    option.inframe(2) = option.max * 2;
end
if isamir(x,'miraudio') && not(isframed(x))
    x = mirframe(x,option.inframe(1),'s',option.inframe(2),'Hz');
end
s = mirspectrum(x,'Power','Terhardt',option.band,'dB','Mask');
type = 'mirspectrum';


function f = main(x,option,postoption)
d = get(x,'Data');
fp = get(x,'FramePos');
fl = option.frame.length.val;
fh = option.frame.hop.val;
if ~fl
    f = mirspectrum(x,'AlongBands','Max',option.max,...
                      'Window',0,'NormalLength',...
                      'Resonance','Fluctuation','MinRes',option.mr);
else
    vb = mirverbose;
    mirverbose(0);
                
    df = cell(1,length(d));
    fp2 = cell(1,length(d));
    p2 = cell(1,length(d));
    for i = 1:length(d)
        df{i} = cell(1,length(d{i}));
        fp2{i} = cell(1,length(d{i}));
        for j = 1:length(d{i})
            dur = fp{i}{j}(1,end) - fp{i}{j}(1,1);
            srj = size(d{i}{j},2) / dur;        % Inner frame rate
            flj = round(fl * srj);
                            % Outer frame length in number of inner frames
            fhj = srj / fh; % Outer hop factor in number of inner frames
            n = floor((dur - fl) * fh ) + 1;    % Number of outer frames
            fp2{i}{j} = zeros(2,n);
            for k = 1:n   % For each outer frame, ...
                st = round( (k-1) * fhj) + 1;
                stend = st + flj - 1;
                dk = d{i}{j}(:,st:stend,:);
                fpk = fp{i}{j}(:,st:stend);
                x2 = set(x,'Data',{{dk}},'FramePos',{{fpk}});
                fk = mirspectrum(x2,'AlongBands','Max',10,'Window',0,...
                                    'NormalLength',...
                                    'Resonance','Fluctuation',...
                                    'MinRes',option.mr);
                dfk = mirgetdata(fk);
                if k == 1
                    df{i}{j} = zeros(size(dfk,1),n,size(dfk,3));
                end
                df{i}{j}(:,k,:) = dfk;
                fp2{i}{j}(:,k) = [fpk(1);fpk(end)];
            end
            p = get(fk,'Pos');
            p2{i}{j} = repmat(p{1}{1},[1 n 1]);
        end
    end
    f = set(fk,'Data',df,'FramePos',fp2,'Pos',p2);
    
    mirverbose(vb);
end

if option.sum
    f = mirsummary(f);
end
f = set(f,'Title','Fluctuation');