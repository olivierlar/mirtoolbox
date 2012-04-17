function varargout = mirchromagram(orig,varargin)
%   c = mirchromagram(x) computes the chromagram, or distribution of energy 
%       along pitches, of the audio signal x.
%       (x can be the name of an audio file as well, or a spectrum, ...)
%   Optional argument:
%       c = mirchromagram(...,'Tuning',t): specifies the central frequency
%           (in Hz.) associated to  chroma C.
%               Default value, t = 261.6256 Hz
%       c = mirchromagram(...,'Wrap',w): specifies whether the chromagram is
%           wrapped or not.
%           w = 1: groups all the pitches belonging to same pitch classes
%               (default value)
%           w = 0: pitches are considered as absolute values.
%       c = mirchromagram(...,'Frame',l,h) orders a frame decomposition of window
%           length l (in seconds) and hop factor h, expressed relatively to
%           the window length. For instance h = 1 indicates no overlap.
%           Default values: l = .2 seconds and h = .05
%       c = mirchromagram(...,'Center'): centers the result.
%       c = mirchromagram(...,'Normal',n): performs a n-norm of the
%           resulting chromagram. Toggled off if n = 0
%               Default value: n = Inf (corresponding to a normalization by
%                   the maximum value).
%       c = mirchromagram(...,'Pitch',p): specifies how to label chromas in
%           the figures.
%               p = 1: chromas are labeled using pitch names (default)
%                   alternative syntax: chromagram(...,'Pitch')
%               p = 0: chromas are labeled using MIDI pitch numbers
%       c = mirchromagram(...,'Triangle'): weight the contribution of each
%           frequency with respect to the distance with the actual
%           frequency of the corresponding chroma.
%       c = mirchromagram(...,'Weight',o): specifies the relative radius of
%           the weighting window, with respect to the distance between
%           frequencies of successive chromas.
%           o = 1: each window begins at the centers of the previous one.
%           o = .5: each window begins at the end of the previous one.
%               (default value)
%       mirchromagram(...,'Min',mi) indicates the lowest frequency taken into
%           consideration in the spectrum computation, expressed in Hz.
%           Default value: 100 Hz. (Gomez, 2006)
%       mirchromagram(...,'Max',ma) indicates the highest frequency taken into
%           consideration in the spectrum computation, expressed in Hz.
%           This upper limit is further shifted to a highest value until
%           the frequency range covers an exact multiple of octaves.
%           Default value: 5000 Hz. (Gomez, 2006)
%       mirchromagram(...,'Res',r) indicates the resolution of the
%           chromagram in number of bins per octave.
%               Default value, r = 12.
%
% Gómez, E. (2006). Tonal description of music audio signal. Phd thesis, 
%   Universitat Pompeu Fabra, Barcelona .

        cen.key = 'Center';
        cen.type = 'Boolean';
        cen.default = 0;
    option.cen = cen;
    
        nor.key = {'Normal','Norm'};
        nor.type = 'Integer';
        nor.default = Inf;
    option.nor = nor;
    
        wth.key = 'Weight';
        wth.type = 'Integer';
        wth.default = .5;
    option.wth = wth;
    
        tri.key = 'Triangle';
        tri.type = 'Boolean';
        tri.default = 0;
    option.tri = tri;
    
        wrp.key = 'Wrap';
        wrp.type = 'Boolean';
        wrp.default = 1;
    option.wrp = wrp;
    
        plabel.key = 'Pitch';
        plabel.type = 'Boolean';
        plabel.default = 1;
    option.plabel = plabel;
    
        thr.key = {'Threshold','dB'};
        thr.type = 'Integer';
        thr.default = 20;
    option.thr = thr;
    
        min.key = 'Min';
        min.type = 'Integer';
        min.default = 100;
    option.min = min;
    
        max.key = 'Max';
        max.type = 'Integer';
        max.default = 5000;
    option.max = max;

        res.key = 'Res';
        res.type = 'Integer';
        res.default = 12;
    option.res = res;

        origin.key = 'Tuning';
        origin.type = 'Integer';
        origin.default = 261.6256;
    option.origin = origin;
    
        transp.key = 'Transpose';
        transp.type = 'Integer';
        transp.default = 0;
    option.transp = transp;
   
specif.option = option;
specif.defaultframelength = .2;
specif.defaultframehop = .05;

varargout = mirfunction(@mirchromagram,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if isamir(x,'mirtemporal') || isamir(x,'mirspectrum')
    freqmin = option.min;
    freqmax = freqmin*2;
    while freqmax < option.max
        freqmax = freqmax*2;
    end
    %freqres = freqmin*(2.^(1/option.res)-1);
        % Minimal frequency resolution should correspond to frequency range
        %   between the first two bins of the chromagram 
        
    x = mirspectrum(x,'dB',option.thr,'Min',freqmin,'Max',freqmax,...
                      'NormalInput','MinRes',option.res,'OctaveRatio',.85);
                  %freqres*.5,...
                  %    'WarningRes',freqres);
end
type = 'mirchromagram';


function c = main(orig,option,postoption)
if iscell(orig)
    orig = orig{1};
end
if option.res == 12
    chromascale = {'C','C#','D','D#','E','F','F#','G','G#','A','A#','B'};
else
    chromascale = 1:option.res;
    option.plabel = 0;
end
if isa(orig,'mirchromagram')
    c = modif(orig,option,chromascale);
else
    c.plabel = 1;
    c.wrap = 0;
    c.chromaclass = {};
    c.chromafreq = {};
    c.register = {};
    c = class(c,'mirchromagram',mirdata(orig));
    c = purgedata(c);
    c = set(c,'Title','Chromagram','Ord','magnitude','Interpolable',0);
    if option.wrp
        c = set(c,'Abs','chroma class');
    else
        c = set(c,'Abs','chroma');
    end
    m = get(orig,'Magnitude');
    f = get(orig,'Frequency');
    %disp('Computing chromagram...')
    fs = get(orig,'Sampling');
    n = cell(1,length(m));  % The final structured list of magnitudes.
    cc = cell(1,length(m));  % The final structured list of chroma classes.
    o = cell(1,length(m));  % The final structured list of octave registers.
    p = cell(1,length(m));  % The final structured list of chromas.
    cf = cell(1,length(m));  % The final structured list of central frequencies related to chromas.
    for i = 1:length(m)
        mi = m{i};
        fi = f{i};
        if not(iscell(mi))
            mi = {mi};
            fi = {fi};
        end
        ni = cell(1,length(mi));    % The list of magnitudes.
        ci = cell(1,length(mi));    % The list of chroma classes.
        oi = cell(1,length(mi));    % The list of octave registers.
        pi = cell(1,length(mi));    % The list of absolute chromas.
        cfi = cell(1,length(mi));    % The central frequency of each chroma.
        for j = 1:length(mi)
            mj = mi{j};
            fj = fi{j};
                        
            % Let's remove the frequencies exceeding the last whole octave.
            minfj = min(min(min(fj)));
            maxfj = max(max(max(fj)));
            maxfj = minfj*2^(floor(log2(maxfj/minfj)));
            fz = find(fj(:,1,1,1) > maxfj);
            mj(fz,:,:,:) = [];      
            fj(fz,:,:,:) = [];
            
            [s1 s2 s3] = size(mj);
            
            cj = freq2chro(fj,option.res,option.origin);
            if not(ismember(min(cj)+1,cj))
                warning('WARNING IN MIRCHROMAGRAM: Frequency resolution of the spectrum is too low.');    
                display('The conversion of low frequencies into chromas may be incorrect.');
            end
            ccj = min(min(min(cj))):max(max(max(cj)));
            sc = length(ccj);   % The size of range of absolute chromas.
            mat = zeros(s1,sc);
            fc = chro2freq(ccj,option.res,option.origin);   % The absolute chromas in Hz.
            fl = chro2freq(ccj-1,option.res,option.origin); % Each previous chromas in Hz.
            fr = chro2freq(ccj+1,option.res,option.origin); % Each related next chromas in Hz.
            for k = 1:sc
                rad = find(and(fj(:,1) > fc(k)-option.wth*(fc(k)-fl(k)),...
                               fj(:,1) < fc(k)-option.wth*(fc(k)-fr(k))));
                if option.tri
                    dist = fc(k) - fj(:,1,1,1);
                    rad1 = dist/(fc(k) - fl(k))/option.wth;
                    rad2 = dist/(fc(k) - fr(k))/option.wth;
                    ndist = max(rad1,rad2);
                    mat(:,k) = max(min(1-ndist,1),0)/length(rad);
                else
                    mat(rad,k) = ones(length(rad),1)/length(rad);
                end
                if k ==1 || k == sc
                    mat(:,k) = mat(:,k)/2;
                end
            end
            nj = zeros(sc,s2,s3);
            for k = 1:s2
                for l = 1:s3
                    nj(:,k,l) = (mj(:,k,l)'*mat)';
                end
            end
            cj = mod(ccj',option.res);
            oi{j} = floor(ccj/option.res)+4;
            if option.plabel
                pj = strcat(chromascale(cj+1)',num2str(oi{j}'));
            else
                pj = ccj'+60;
            end
            ci{j} = repmat(cj,[1,s2,s3]);
            pi{j} = repmat(pj,[1,s2,s3]);
            ni{j} = nj;
            cfi{j} = fc;
        end
        n{i} = ni;
        cc{i} = ci;
        o{i} = oi;
        p{i} = pi;
        cf{i} = cfi;
    end
    c = set(c,'Magnitude',n,'Chroma',p,'ChromaClass',cc,...
              'ChromaFreq',cf,'Register',o);
    c = modif(c,option,chromascale);
    c = {c orig};
end

   
function c = modif(c,option,chromascale)
if option.plabel
    c = set(c,'PitchLabel',1);
end            
if option.cen || option.nor || option.wrp || option.transp
    n = get(c,'Magnitude');
    p = get(c,'Chroma');
    cl = get(c,'ChromaClass');
    fp = get(c,'FramePos');
    n2 = cell(1,length(n));
    p2 = cell(1,length(n));
    if option.transp
        transp = mod(option.transp,12);
    end
    wrp = option.wrp && not(get(c,'Wrap'));
    for i = 1:length(n)
        ni = n{i};
        pi = p{i};
        cli = cl{i};
        if not(iscell(ni))
            ni = {ni};
            pi = {pi};
            cli = {cli};
        end
        if wrp
            c = set(c,'Wrap',option.wrp);
        end
        n2i = cell(1,length(ni));
        p2i = cell(1,length(ni));
        for j = 1:length(ni)
            nj = ni{j};
            pj = pi{j};
            clj = cli{j};
            if wrp
                n2j = zeros(option.res,size(nj,2),size(nj,3));
                for k = 1:size(pj,1)
                    n2j(clj(k)+1,:,:) = n2j(clj(k)+1,:,:) + nj(k,:,:);  % squared sum (parameter)
                end
                p2i{j} = chromascale';
            else
                n2j = nj;
                p2i{j} = pi{j};
            end
            if option.cen
                n2j = n2j - repmat(mean(n2j),[size(n2j,1),1,1]);
            end
            if option.nor
                n2j = n2j ./ repmat(vectnorm(n2j,option.nor) + ...
                    repmat(1e-6,[1,size(n2j,2),size(n2j,3)] )...
                    ,[size(n2j,1),1,1]);
            end
            if option.transp
                n2j = [n2j(13-transp:end,:,:);n2j(1:12-transp,:,:)];
            end
            n2i{j} = n2j;
        end
        n2{i} = n2i;
        p2{i} = p2i;
    end
    c = set(c,'Magnitude',n2,'Chroma',p2,'FramePos',fp);
end


function c = freq2chro(f,res,origin)
c = round(res*log2(f/origin));


function f = chro2freq(c,res,origin)
f = 2.^(c/res)*origin;


function y = vectnorm(x,p)
if isinf(p)
    y = max(x);
else
    y = sum(abs(x).^p).^(1/p);
end