function varargout = mirenvelope(orig,varargin)
%   e = mirenvelope(x) extracts the envelope of x, showing the global shape
%       of the waveform.
%   mirenvelope(...,m) specifies envelope extraction method.
%       Possible values:
%           m = 'Filter' uses a low-pass filtering. (Default strategy)
%           m = 'Spectro' uses a spectrogram.
%
%   Options related to the 'Filter' method:
%       mirenvelope(...,'Hilbert'): performs a preliminary Hilbert
%           transform.
%       mirenvelope(...,'PreDecim',N) downsamples by a factor N>1, where 
%           N is an integer, before the low-pass filtering (Klapuri, 1999).
%           Default value: N = 1.
%       mirenvelope(...,'Filtertype',f) specifies the filter type.
%           Possible values are:
%               f = 'IIR': filter with one autoregressive coefficient 
%                   (default)
%               f = 'HalfHann': half-Hanning (raised cosine) filter 
%                   (Scheirer, 1998)
%           Option related to the 'IIR' option:
%           mirenvelope(...,'Tau',t): time constant of low-pass filter in 
%               seconds.
%               Default value: t = 0.02 s.
%       mirenvelope(...,'PostDecim',N) downsamples  by a factor N>1, where 
%           N is an integer, after the low-pass filtering.
%           Default value: N = 16 if 'PreDecim' is not used, else N = 1.
%       mirenvelope(...,'Trim'): trims the initial ascending phase of the
%           curves related to the transitory state.
%
%   Options related to the 'Spectro' method:
%       mirenvelope(...,b) specifies whether the frequency range is further
%           decomposed into bands. Possible values:
%               b = 'Freq': no band decomposition (default value)
%               b = 'Mel': Mel-band decomposition 
%               b = 'Bark': Bark-band decomposition 
%               b = 'Cents': decompositions into cents 
%       mirenvelope(...,'Frame',...) specifies the frame configuration.
%           Default value: length: .1 s, hop factor: 10 %.
%       mirenvelope(...,'UpSample',N) upsamples by a factor N>1, where 
%           N is an integer.
%           Default value if 'UpSample' called: N = 2
%       mirenvelope(...,'Complex') toggles on the 'Complex' method for the
%           spectral flux computation.
%
%   Other available for all methods:
%       mirenvelope(...,'Sampling',r): resamples to rate r (in Hz).
%           'Down' and 'Sampling' options cannot therefore be combined.
%       mirenvelope(...,'Halfwave'): performs a half-wave rectification.
%       mirenvelope(...,'Center'): centers the extracted envelope.
%       mirenvelope(...,'HalfwaveCenter'): performs a half-wave
%           rectification on the centered envelope.
%       mirenvelope(...,'Log',mu): computes the logarithm of the
%           envelope, before the eventual differentiation, using a mu-law
%           compression (Klapuri, 2006).
%               Default value for mu: 100
%       mirenvelope(...,'Power'): computes the power (square) of the
%           envelope.
%       mirenvelope(...,'Diff'): computes the differentation of the
%           envelope, i.e., the differences between successive samples.
%       mirenvelope(...,'HalfwaveDiff'): performs a half-wave
%           rectification on the differentiated envelope.
%       mirenvelope(...,'Normal'): normalizes the values of the envelope by
%           fixing the maximum value to 1.
%       mirenvelope(...,'Lambda',l): sums the half-wave rectified envelope
%           with the non-differentiated envelope, using the respective
%           weight 0<l<1 and (1-l). (Klapuri et al., 2006)
%       mirenvelope(...,'Smooth',o): smooths the envelope using a moving
%           average of order o.
%           Default value when the option is toggled on: o=30
%       mirenvelope(...,'Gauss',o): smooths the envelope using a gaussian
%           of standard deviation o samples.
%           Default value when the option is toggled on: o=30
%       mirenvelope(...,'Klapuri06'): follows the model proposed in
%           (Klapuri et al., 2006). 

        method.type = 'String';
        method.choice = {'Filter','Spectro'};
        method.default = 'Filter';
    option.method = method;
    
%% options related to 'Filter':

        hilb.key = 'Hilbert';
        hilb.type = 'Boolean';
        hilb.default = 0;
    option.hilb = hilb;
 
        decim.key = {'Decim','PreDecim'};
        decim.type = 'Integer';
        decim.default = 0;
    option.decim = decim;
    
        filter.key = 'FilterType';
        filter.type = 'String';
        filter.choice = {'IIR','HalfHann',0};
        if isamir(orig,'mirenvelope')
            filter.default = 0; % no more envelope extraction, already done
        else
            filter.default = 'IIR';
        end
    option.filter = filter;

            %% options related to 'IIR': 
            tau.key = 'Tau';
            tau.type = 'Integer';
            tau.default = .02;
    option.tau = tau;
    
        zp.key = 'ZeroPhase'; % internal use: for manual filtfilt
        zp.type = 'Boolean';
        if isamir(orig,'mirenvelope')
            zp.default = 0;
        else
            zp.default = NaN;
        end
    option.zp = zp;

        ds.key = {'Down','PostDecim'};
        ds.type = 'Integer';
        if isamir(orig,'mirenvelope')
            ds.default = 1;
        else
            ds.default = NaN; % 0 if 'PreDecim' is used, else 16
        end
        ds.when = 'After';
        ds.chunkcombine = 'During';
    option.ds = ds;

        trim.key = 'Trim';
        trim.type = 'Boolean';
        trim.default = 0;
        trim.when = 'After';
    option.trim = trim;
     
%% Options related to 'Spectro':

        band.type = 'String';
        band.choice = {'Freq','Mel','Bark','Cents'};
        band.default = 'Freq';
    option.band = band;

        up.key = {'UpSample'};
        up.type = 'Integer';
        up.default = 0;
        up.keydefault = 2;
        up.when = 'After';
    option.up = up;

        complex.key = 'Complex';
        complex.type = 'Boolean';
        complex.default = 0;
        complex.when = 'After';
    option.complex = complex;    
    
%% Options related to all methods:
    
        sampling.key = 'Sampling';
        sampling.type = 'Integer';
        sampling.default = 0;
        sampling.when = 'After';
    option.sampling = sampling;

        hwr.key = 'Halfwave';
        hwr.type = 'Boolean';
        hwr.default = 0;
        hwr.when = 'After';
    option.hwr = hwr;

        c.key = 'Center';
        c.type = 'Boolean';
        c.default = 0;
        c.when = 'After';
    option.c = c;
    
        chwr.key = 'HalfwaveCenter';
        chwr.type = 'Boolean';
        chwr.default = 0;
        chwr.when = 'After';
    option.chwr = chwr;
    
        oplog.key = 'Log';
        oplog.type = 'Integer';
        oplog.default = 0;
        oplog.keydefault = 100;
        oplog.when = 'After';
    option.log = oplog;
    
        oppow.key = 'Power';
        oppow.type = 'Integer';
        oppow.default = 0;
        oppow.when = 'After';
    option.power = oppow;

        diff.key = 'Diff';
        diff.type = 'Integer';
        diff.default = 0;
        diff.keydefault = 1;
        diff.when = 'After';
    option.diff = diff;
    
        diffhwr.key = 'HalfwaveDiff';
        diffhwr.type = 'Integer';
        diffhwr.default = 0;
        diffhwr.keydefault = 1;
        diffhwr.when = 'After';
    option.diffhwr = diffhwr;

        lambda.key = 'Lambda';
        lambda.type = 'Integer';
        lambda.default = 1;
        lambda.when = 'After';
    option.lambda = lambda;

        aver.key = 'Smooth';
        aver.type = 'Integer';
        aver.default = 0;
        aver.keydefault = 30;
        aver.when = 'After';
    option.aver = aver;
        
        gauss.key = 'Gauss';
        gauss.type = 'Integer';
        gauss.default = 0;
        gauss.keydefault = 30;
        gauss.when = 'After';
    option.gauss = gauss;

        norm.key = 'Normal';
        norm.type = 'Boolean';
        norm.default = 0;
        norm.when = 'After';
    option.norm = norm;

        presel.type = 'String';
        presel.choice = {'Klapuri06'};
        presel.default = 0;
    option.presel = presel;    
    
        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [.1 .1];
    option.frame = frame;
        
specif.option = option;

specif.eachchunk = 'Normal';
specif.combinechunk = 'Concat';
specif.extensive = 1;

varargout = mirfunction(@mirenvelope,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
type = 'mirenvelope';
if isamir(x,'mirscalar')
    return
end
if ischar(option.presel) && strcmpi(option.presel,'Klapuri06')
    option.method = 'Spectro';
end
if not(isamir(x,'mirenvelope'))
    if strcmpi(option.method,'Filter')
        if isnan(option.zp)
            if strcmpi(option.filter,'IIR')
                option.zp = 1;
            else
                option.zp = 0;
            end
        end
        if option.zp == 1
            x = mirenvelope(x,'ZeroPhase',2,'Down',1,...
                              'Tau',option.tau,'PreDecim',option.decim);
        end
    elseif strcmpi(option.method,'Spectro')
        x = mirspectrum(x,'Frame',option.frame.length.val,...
                                  option.frame.length.unit,...
                                  option.frame.hop.val,...
                                  option.frame.hop.unit,...
                                  'Window','hanning',...
                                  option.band,'Power');

    end
end


function e = main(orig,option,postoption)
if iscell(orig)
    orig = orig{1};
end
if isamir(orig,'mirscalar')
    d = get(orig,'Data');
    fp = get(orig,'FramePos');
    for i = 1:length(d)
        for j = 1:length(d{i})
            d{i}{j} = reshape(d{i}{j},size(d{i}{j},2),1,size(d{i}{j},3));
            p{i}{j} = mean(fp{i}{j})';
        end
    end
    e.downsampl = 0;
    e.hwr = 0;
    e.diff = 0;
    e.method = 'Spectro';
    e.phase = {{}};
    e = class(e,'mirenvelope',mirtemporal(orig));
    e = set(e,'Title','Envelope','Data',d,'Pos',p,'FramePos',{{}});
    postoption.trim = 0;
    postoption.ds = 0;
    e = post(e,postoption);
    return
end
if isfield(option,'presel') && ischar(option.presel) && ...
        strcmpi(option.presel,'Klapuri06')
    option.method = 'Spectro';
    postoption.up = 2;
    postoption.log = 100;
    postoption.diffhwr = 1;
    postoption.lambda = .8;
end
if isfield(postoption,'ds') && isnan(postoption.ds)
    if option.decim
        postoption.ds = 0;
    else
        postoption.ds = 16;
    end
end
if not(isfield(option,'filter')) || not(ischar(option.filter))
    e = post(orig,postoption);
elseif strcmpi(option.method,'Spectro')
    d = get(orig,'Data');
    fp = get(orig,'FramePos');
    sr = get(orig,'Sampling');
    ch = get(orig,'Channels');
    ph = get(orig,'Phase');
    for h = 1:length(d)
        for i = 1:length(d{h})
            if size(d{h}{i},3)>1 % Already in bands (channels in 3d dim)
                d{h}{i} = permute(sum(d{h}{i}),[2 1 3]);
                ph{h}{i} = permute(ph{h}{i},[2 1 3]);
            else % Simple spectrogram, frequency range sent to 3d dim
                d{h}{i} = permute(d{h}{i},[2 3 1]);
                ph{h}{i} = permute(ph{h}{i},[2 3 1]);
            end
            p{h}{i} = mean(fp{h}{i})';
        end
        if size(fp{h}{i},2)<2
            %error('ERROR IN MIRENVELOPE: The frame decomposition did not succeed. Either the input is of too short duration, or the chunk size is too low.');
            sr{h} = 0;
        else
            sr{h} = 1/(fp{h}{i}(1,2)-fp{h}{i}(1,1));
        end
        ch{h} = (1:size(d{h}{i},3))';
    end
    e.downsampl = 0;
    e.hwr = 0;
    e.diff = 0;
    %e.todelete = {};
    e.method = 'Spectro';
    e.phase = ph;
    e = class(e,'mirenvelope',mirtemporal(orig));
    e = set(e,'Title','Envelope','Data',d,'Pos',p,...
              'Sampling',sr,'Channels',ch);
    postoption.trim = 0;
    postoption.ds = 0;
    e = post(e,postoption);
else
    if isnan(option.zp)
        if strcmpi(option.filter,'IIR')
            option.zp = 1;
        else
            option.zp = 0;
        end
    end
    if option.zp == 1
        option.decim = 0;
    end
    e.downsampl = 1;
    e.hwr = 0;
    e.diff = 0;
    e.method = option.filter;
    e.phase = {};
    e = class(e,'mirenvelope',mirtemporal(orig));
    e = purgedata(e);
    e = set(e,'Title','Envelope');
    sig = get(e,'Data');
    x = get(e,'Pos');
    sr = get(e,'Sampling');
    disp('Extracting envelope...')
    d = cell(1,length(sig));
    [state e] = gettmp(orig,e);
    for k = 1:length(sig)
        if option.decim
            sr{k} = sr{k}/option.decim;
        end
        if strcmpi(option.filter,'IIR')
            a2 = exp(-1/(option.tau*sr{k})); % filter coefficient 
            a = [1 -a2];
            b = 1-a2;
        elseif strcmpi(option.filter,'HalfHann')
            a = 1;
            b = hann(sr{k}*.4);
            b = b(ceil(length(b)/2):end);
        end
        d{k} = cell(1,length(sig{k}));
        for i = 1:length(sig{k})
            sigi = sig{k}{i};
            if option.zp == 2
                sigi = flipdim(sigi,1);
            end
            if option.hilb
                try
                    for h = 1:size(sigi,2)
                        for j = 1:size(sigi,3)
                            sigi(:,h,j) = hilbert(sigi(:,h,j));
                        end
                    end
                catch
                    disp('Signal Processing Toolbox does not seem to be installed. No Hilbert transform.');
                end 
            end
            sigi = abs(sigi);
            
            if option.decim
                dsigi = zeros(ceil(size(sigi,1)/option.decim),...
                              size(sigi,2),size(sigi,3));
                for f = 1:size(sigi,2)
                    for c = 1:size(sigi,3)
                        dsigi(:,f,c) = decimate(sigi(:,f,c),option.decim);
                    end
                end
                sigi = dsigi;
                clear dsigi
                x{k}{i} = x{k}{i}(1:option.decim:end,:,:);
            end
            
            % tmp = filtfilt(1-a,[1 -a],sigi); % zero-phase IIR filter for smoothing the envelope

            % Manual filtfilt
            emptystate = isempty(state);
            tmp = zeros(size(sigi));
            for c = 1:size(sigi,3)
                if emptystate
                    [tmp(:,:,c) state(:,c,1)] = filter(b,a,sigi(:,:,c));
                else
                    [tmp(:,:,c) state(:,c,1)] = filter(b,a,sigi(:,:,c),...
                                                        state(:,c,1));
                end
            end
            
            tmp = max(tmp,0); % For security reason...
            if option.zp == 2
                tmp = flipdim(tmp,1);
            end
            d{k}{i} = tmp;
            %td{k} = round(option.tau*sr{k}*1.5); 
        end
    end
    e = set(e,'Data',d,'Pos',x,'Sampling',sr); %,'ToDelete',td
    e = settmp(e,state);
    if not(option.zp == 2)
        e = post(e,postoption);
    end
end
if isfield(option,'presel') && ischar(option.presel) && ...
        strcmpi(option.presel,'Klapuri06')
    e = mirsum(e,'Adjacent',10);
end

    
function e = post(e,postoption)
if isempty(postoption)
    return
end
if isfield(postoption,'lambda') && not(postoption.lambda)
    postoption.lambda = 1;
end
d = get(e,'Data');
tp = get(e,'Time');
sr = get(e,'Sampling');
ds = get(e,'DownSampling');
ph = get(e,'Phase');
for k = 1:length(d)
    if isfield(postoption,'sampling')
        if postoption.sampling
            newsr = postoption.sampling;
        elseif isfield(postoption,'ds') && postoption.ds>1
            newsr = sr{k}/postoption.ds;
        else
            newsr = sr{k};
        end
    end
    if isfield(postoption,'up') && postoption.up
        [z,p,gain] = butter(6,10/newsr/postoption.up*2,'low');
        [sos,g] = zp2sos(z,p,gain);
        Hd = dfilt.df2tsos(sos,g);
    end
    for i = 1:length(d{k})
        if isfield(postoption,'sampling') && postoption.sampling
            if and(sr{k}, not(sr{k} == postoption.sampling))
                dk = d{k}{i};
                for j = 1:size(dk,3)
                    if not(sr{k} == round(sr{k}))
                        mirerror('mirenvelope','The ''Sampling'' postoption cannot be used after using the ''Down'' postoption.');
                    end
                    rk(:,:,j) = resample(dk(:,:,j),postoption.sampling,sr{k});
                end
                d{k}{i} = rk;
                tp{k}{i} = repmat((0:size(d{k}{i},1)-1)',...
                                    [1 1 size(tp{k}{i},3)])...
                            /postoption.sampling + tp{k}{i}(1,:,:);
                if not(iscell(ds))
                    ds = cell(length(d));
                end
                ds{k} = round(sr{k}/postoption.sampling);
            end
        elseif isfield(postoption,'ds') && postoption.ds>1
            if not(postoption.ds == round(postoption.ds))
                mirerror('mirenvelope','The ''Down'' sampling rate should be an integer.');
            end
            ds = postoption.ds;
            tp{k}{i} = tp{k}{i}(1:ds:end,:,:); % Downsampling...
            d{k}{i} = d{k}{i}(1:ds:end,:,:);
        end
        if isfield(postoption,'sampling')
            if not(strcmpi(e.method,'Spectro')) && postoption.trim 
                tdk = round(newsr*.1); 
                d{k}{i}(1:tdk,:,:) = repmat(d{k}{i}(tdk,:,:),[tdk,1,1]); 
                d{k}{i}(end-tdk+1:end,:,:) = repmat(d{k}{i}(end-tdk,:,:),[tdk,1,1]);
            end
            if postoption.log
                dki = max(0,d{k}{i});
                mu = postoption.log;
                dki = log(1+mu*dki)/log(1+mu);
                dki(~isfinite(d{k}{i})) = NaN;
                d{k}{i} = dki;
            end
            if postoption.power
                d{k}{i} = d{k}{i}.^2;
            end
            if postoption.up
                dki = zeros(size(d{k}{i},1).*postoption.up,...
                            size(d{k}{i},2),size(d{k}{i},3));
                dki(1:postoption.up:end,:,:) = d{k}{i};
                dki = filter(Hd,[dki;...
                                zeros(6,size(d{k}{i},2),size(d{k}{i},3))]);
                d{k}{i} = dki(1+ceil(6/2):end-floor(6/2),:,:);
                tki = zeros(size(tp{k}{i},1).*postoption.up,...
                            size(tp{k}{i},2),...
                            size(tp{k}{i},3));
                dt = repmat((tp{k}{i}(2)-tp{k}{i}(1))...
                                /postoption.up,...
                            [size(tp{k}{i},1),1,1]);
                for j = 1:postoption.up
                    tki(j:postoption.up:end,:,:) = tp{k}{i}+dt*(j-1);
                end
                tp{k}{i} = tki;
                newsr = sr{k}*postoption.up;
            end
            if (postoption.diffhwr || postoption.diff) && ...
                    not(get(e,'Diff'))
                tp{k}{i} = tp{k}{i}(1:end-1,:,:);
                order = max(postoption.diffhwr,postoption.diff);
                if postoption.complex
                    dph = diff(ph{k}{i},2);
                    dph = dph/(2*pi);% - round(dph/(2*pi));
                    ddki = sqrt(d{k}{i}(3:end,:,:).^2 + d{k}{i}(2:end-1,:,:).^2 ...
                                              - 2.*d{k}{i}(3:end,:,:)...
                                                 .*d{k}{i}(2:end-1,:,:)...
                                                 .*cos(dph));
                    d{k}{i} = d{k}{i}(2:end,:,:); 
                    tp{k}{i} = tp{k}{i}(2:end,:,:);
                elseif order == 1
                    ddki = diff(d{k}{i});
                else
                    b = firls(order,[0 0.9],[0 0.9*pi],'differentiator');
                    ddki = filter(b,1,...
                        [repmat(d{k}{i}(1,:,:),[order,1,1]);...
                         d{k}{i};...
                         repmat(d{k}{i}(end,:,:),[order,1,1])]);
                    ddki = ddki(order+1:end-order-1,:,:);
                end
                if postoption.diffhwr
                    ddki = hwr(ddki);
                end
                d{k}{i} = (1-postoption.lambda)*d{k}{i}(1:end-1,:,:)...
                            + postoption.lambda*sr{k}/10*ddki;
            end
            if postoption.aver
                y = filter(ones(1,postoption.aver),1,...
                            [d{k}{i};zeros(postoption.aver,...
                                           size(d{k}{i},2),...
                                           size(d{k}{i},3))]);
                d{k}{i} = y(1+ceil(postoption.aver/2):...
                             end-floor(postoption.aver/2),:,:);
            end
            if postoption.gauss
                sigma = postoption.gauss;
                gauss = 1/sigma/2/pi...
                        *exp(- (-4*sigma:4*sigma).^2 /2/sigma^2);
                y = filter(gauss,1,[d{k}{i};zeros(4*sigma,1,size(d{k}{i},3))]);
                y = y(4*sigma:end,:,:);
                d{k}{i} = y(1:size(d{k}{i},1),:,:);
            end
            if postoption.chwr
                d{k}{i} = center(d{k}{i});
                d{k}{i} = hwr(d{k}{i});
            end
            if postoption.hwr
                d{k}{i} = hwr(d{k}{i});
            end
            if postoption.c
                d{k}{i} = center(d{k}{i});
            end    
            if postoption.norm
                d{k}{i} = d{k}{i}./repmat(max(abs(d{k}{i})),...
                                         [size(d{k}{i},1),1,1]);
            end        
        end
    end
    if isfield(postoption,'sampling')
        sr{k} = newsr;
    end
end
if isfield(postoption,'ds') && postoption.ds>1
    e = set(e,'DownSampling',postoption.ds,'Sampling',sr);
elseif isfield(postoption,'sampling') && postoption.sampling
    e = set(e,'DownSampling',ds,'Sampling',sr);
elseif isfield(postoption,'up') && postoption.up
    e = set(e,'Sampling',sr);
end
if isfield(postoption,'sampling')
    if postoption.hwr
        e = set(e,'Halfwave',1);
    end
    if postoption.diff
        e = set(e,'Diff',1,'Halfwave',0,'Title','Differentiated envelope');
    end
    if postoption.diffhwr
        e = set(e,'Diff',1,'Halfwave',1,'Centered',0);
    end
    if postoption.c
        e = set(e,'Centered',1);
    end
    if postoption.chwr
        e = set(e,'Halfwave',1,'Centered',1);
    end
end
e = set(e,'Data',d,'Time',tp); 
%if isfield(postoption,'frame') && isstruct(postoption.frame)
%    e = mirframenow(e,postoption);
%end