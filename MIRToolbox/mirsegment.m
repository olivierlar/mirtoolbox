function [f,p,m,fe] = mirsegment(x,varargin)
%   f = mirsegment(a) segments an audio signal. It can also be the name of an
%       audio file or 'Folder', for the analysis of the audio files in the
%       current folder. The segmentation of audio signal already decomposed
%       into frames is not available for the moment.
%   f = mirsegment(...,strategy) segments following alternative strategies:
%       'Novelty': using a self-similarity matrix (Foote & Cooper, 2003).
%           (by default)
%       'HCDF': using the Harmonic Change Detection Function (Harte &
%           Sandler, 2006)
%       The options available for the chosen strategies can be specified
%           directly as options of the segment function.
%           Example: mirsegment(a,'Novelty','KernelSize',10)
%   f = mirsegment(...,feature) base the segmentation strategy on a specific
%       feature.
%       'Spectrum': from the spectrum of the audio signal x.
%           (by default)
%       'MFCC': from the MFCC of the audio signal x.
%       'Keystrength': from the key strength profile of the audio signal x.
%       The option related to this feature extraction can be specified too.
%       Example: mirsegment(...,'Spectrum','Window','bartlett')
%                mirsegment(...,'MFCC','Rank',1:10)
%                mirsegment(...,'Keystrength','Weight',.5)
%   These feature need to be frame-based, in order to appreciate their
%       temporal evolution. Therefore, the audio signal x is first
%       decomposed into frames. This decomposition can be controled using
%       the 'Frame' keyword.  
%
%   f = mirsegment(a,s) segments a using the results of a segmentation
%       analysis s. s can be the peaks detected on an analysis of the
%       audio for instance.
%
%   f = mirsegment(a,v) where v is an array of numbers, segments a using
%       the temporal positions specified in v (in s.)
%
%   [f,p] = mirsegment(...) also displays the analysis produced by the chosen
%       strategy.
%           For 'Novelty', p is the novelty curve.
%           For 'HCDF', p is the Harmonic Change Detection Function.
%   [f,p,m] = mirsegment(...) also displays the preliminary analysis
%       undertaken in the chosen strategy.
%           For 'Novelty', m is the similarity matrix.
%           For 'HCDF', m is the tonal centroid.
%   [f,p,m,fe] = mirsegment(...) also displays the temporal evolution of the
%       feature used for the analysis.
%
%   Foote, J. & Cooper, M. (2003). Media Segmentation using Self-Similarity
%       Decomposition,. In Proc. SPIE Storage and Retrieval for Multimedia
%       Databases, Vol. 5021, pp. 167-75.
%   Harte, C. A. & Sandler, M. B. (2006). Detecting harmonic change in
%       musical audio, in Proceedings of Audio and Music Computing for 
%       Multimedia Workshop, Santa Barbara, CA.
 
        ana.type = 'String';
        ana.choice = {'Spectrum','Keystrength','Pitch'};
        ana.default = 0;
    option.ana = ana;
    
        band.choice = {'Mel','Bark','Freq'};
        band.type = 'String';
        band.default = 'Freq';
    option.band = band;

        mfc.key = {'Rank','MFCC'};
        mfc.type = 'Integer';
        mfc.default = 0;
        mfc.keydefault = 1:13;
    option.mfc = mfc;

        strat.choice = {'Novelty','HCDF'};
        strat.default = 'Novelty';
        strat.position = 2;
    option.strat = strat;

        K.key = 'KernelSize';
        K.type = 'Integer';
        K.default = 128;
    option.K = K;
    
        distance.key = 'Distance';
        distance.type = 'String';
        distance.default = 'cosine';
    option.distance = distance;

        measure.key = {'Measure','Similarity'};
        measure.type = 'String';
        measure.default = 'exponential';
    option.measure = measure;

        tot.key = 'Total';
        tot.type = 'Integer';
        tot.default = Inf;
    option.tot = tot;

        cthr.key = 'Contrast';
        cthr.type = 'Integer';
        cthr.default = .1;
    option.cthr = cthr;
    
        mi.key = 'Min';
        mi.type = 'Integer';
        mi.default = 0;
    option.mi = mi;
     
        ma.key = 'Max';
        ma.type = 'Integer';
        ma.default = 0;
    option.ma = ma;
    
        norm.key = 'Normal';
        norm.type = 'Boolean';
        norm.default = 0;
    option.norm = norm;

        win.key = 'Window';
        win.type = 'String';
        win.default = 'hamming';
    option.win = win;

        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [3 .1];
    option.frame = frame;
%defaultframelength = 0.05;   %.5 for keystrength, .743 for HCDF
%defaultframehop = 1;         %.2 for keystrength, .01 for pitch, 1/8 for hcdf
   
specif.option = option;


p = {};
m = {};
fe = {};

if isa(x,'mirdesign')
    if not(get(x,'Eval'))
        % During bottom-up construction of the general design

        [unused option] = miroptions(@mirframe,x,specif,varargin);
        type = get(x,'Type');
        f = mirdesign(@mirsegment,x,option,{},struct,type);
        
        sg = get(x,'Segment');
        if not(isempty(sg))
            f = set(f,'Segment',sg);
        else
            f = set(f,'Segment',option.strat);
        end
        
    else
        % During top-down evaluation initiation
        
        e = evaleach(x);
        if iscell(e)
            e = e{1};
        end
        f = mirsegment(e,varargin{:});
        p = x;
    end
elseif isa(x,'mirdata')
    [unused option] = miroptions(@mirframe,x,specif,varargin);
    if ischar(option.strat)
        if strcmpi(option.strat,'Novelty') || strcmpi(option.strat,'HCDF')
            dx = get(x,'Data');
            if size(dx{1},2) > 1
                %fr = x;
                error('ERROR IN MIRSEGMENT: The segmentation of audio signal already decomposed into frames is not available for the moment.');
            elseif not(strcmpi(option.ana,'Pitch'))
                if not(option.frame.length.val)
                    if strcmpi(option.ana,'Keystrength')
                        option.frame.length.val = .5;
                        option.frame.hop.val = .2;
                    elseif strcmpi(option.strat,'HCDF')
                        option.frame.length.val = .743;
                        option.frame.hop.val = 1/8;
                    elseif strcmpi(option.ana,'Pitch')
                        option.frame.length.val = .05;
                        option.frame.hop.val = .01;
                    else
                        option.frame.length.val = .05;
                        option.frame.hop.val = 1;
                    end
                end
                fr = mirframenow(x,option);
            end
            if not(isequal(option.mfc,0))
                fe = mirmfcc(fr,'Rank',option.mfc);
            elseif strcmpi(option.ana,'Spectrum')
                fe = mirspectrum(fr,'Min',option.mi,'Max',option.ma,...
                                    'Normal',option.norm,option.band,...
                                    'Window',option.win);
            elseif strcmpi(option.ana,'Keystrength')
                    fe = mirkeystrength(fr);
            elseif strcmpi(option.ana,'Pitch')
                [unused,fe] = mirpitch(x,'Frame');
            else
                fe = fr;
            end
            if strcmpi(option.strat,'Novelty')
                [n,m] = mirnovelty(fe,'Distance',option.distance,...
                                      'Measure',option.measure,...
                                      'KernelSize',option.K);
                p = mirpeaks(n,'Total',option.tot,...
                                'Contrast',option.cthr,...
                                'Chrono','NoBegin','NoEnd');
            elseif strcmpi(option.strat,'HCDF')
                [df m fe] = mirhcdf(fe);
                p = mirpeaks(df);
            end
        else
            error('ERROR IN MIRSEGMENT: Syntax error. See help mirsegment.');
        end
        f = mirsegment(x,p);
    else
        dx = get(x,'Data');
        dt = get(x,'Time');

        if isa(option.strat,'mirscalar')
            ds = get(option.strat,'PeakPos');
            fp = get(option.strat,'FramePos');
        elseif isa(option.strat,'mirdata')
            ds = get(option.strat,'AttackPos');
            if isempty(ds)
                ds = get(option.strat,'PeakPos');
            end
            xx = get(option.strat,'Pos');
        else
            ds = option.strat;
            fp = cell(1,length(dx));
        end
        st = cell(1,length(dx));
        sx = cell(1,length(dx));
        cl = cell(1,length(dx));
        for k = 1:length(dx)
            dxk = dx{k}{1}; % values in kth audio file
            dtk = dt{k}{1}; % time positions in kth audio file
            if isa(option.strat,'mirdata')
                dsk = ds{k}{1}; % segmentation times in kth audio file
            else
                dsk = {ds};
            end
            fsk = [];   % the structured array of segmentation times 
                         % needs to be flatten
            for j = 1:length(dsk)
                if isa(option.strat,'mirdata')
                    dsj = dsk{j}; % segmentation times in jth segment
                else
                    dsj = ds;
                end
                if not(iscell(dsj))
                    dsj = {dsj};
                end
                for m = 1:length(dsj)
                    % segmentation times in mth bank channel
                    if isa(option.strat,'mirscalar')
                        dsm = fp{k}{m}(1,dsj{m});
                    elseif isa(option.strat,'mirdata')
                        dsm = xx{k}{m}(dsj{m});
                    else
                        dsm = dsj{m};
                    end
                    dsm(find(dsm <= dtk(1))) = [];
                    dsm(find(dsm >= dtk(end))) = [];
                    % It is presupposed here that the segmentations times
                    % for a given channel are not decomposed per frames,
                    % because the segmentation of the frame decomposition
                    % is something that does not seem very clear.
                    % Practically, the peak picking for instance is based 
                    % therefore on a frame analysis (such as novelty), and
                    % segmentation are inferred between these frames...
                    fsk = [fsk dsm];
                end
            end

            fsk = sort(fsk); % Here is the chronological ordering

            if isempty(fsk)
                ffsk = {[0;dtk(end)]};
            else
                ffsk = cell(1,length(fsk)+1);
                ffsk{1} = [dtk(1);fsk(1)];
                for h = 1:length(fsk)-1
                    ffsk{h+1} = [fsk(h);fsk(h+1)];
                end
                ffsk{end} = [fsk(end);dtk(end)];
            end
            n = length(ffsk);

            crd = zeros(1,n+1); % the sample positions of the
                                % segmentations in the channel
            crd0 = 0;
            for i = 1:n
                crd0 = crd0 + find(dtk(crd0+1:end)>=ffsk{i}(1),1);
                crd(i) = crd0;
            end
            crd(n+1) = size(dxk,1)+1;

            sxk = cell(1,n); % each cell contains a segment
            stk = cell(1,n); % each cell contains
                             % the corresponding time positions

            for i = 1:n
                sxk{i} = dxk(crd(i):crd(i+1)-1,1,:);
                stk{i} = dtk(crd(i):crd(i+1)-1);
            end
            sx{k} = sxk;
            st{k} = stk;
            fp{k} = ffsk;
            cl{k} = 1:n;
        end
        f = set(x,'Data',sx,'Time',st,'FramePos',fp,'Clusters',cl);
        p = strat;
        m = {};
        fe = {};
    end
else
    f = mirsegment(miraudio(x),varargin{:});
end 