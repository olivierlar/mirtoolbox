function varargout = mirlowenergy(x,varargin)
%   p = mirlowenergy(f) computes the percentage of frames showing a RMS 
%       energy that is lower than a given threshold. 
%   For instance, for a musical excerpt with some very loud frames and 
%       lots of silent frames, we would get a high low-energy rate.
%   Optional argument:
%       mirlowenergy(...,'Threshold',t) expressed as a ratio to the average
%           energy over the frames.
%           Default value: t = 1
%       mirlowenergy(...,'Frame',l,h) specifies the use of frames of
%           length l seconds and a hop rate h.
%           Default values: l = .05 s, h = .5
%       mirlowenergy(...,'Root',0) uses mean square instead of root mean
%           square
%       mirlowenergy(...,'ASR') computes the Average Silence Ratio, which
%           corresponds in fact to mirlowenergy(...,'Root',0,'Threshold',t)
%           where t is fixed here by default to t = .5
%   [p,e] = mirlowenergy(...) also returns the RMS energy curve.

%        nbframes.type = 'Boolean';
%        nbframes.default = 0;
%    option.nbframes = nbframes;
    
        asr.key = 'ASR';
        asr.type = 'Boolean';
        asr.default = 0;
    option.asr = asr;
    
        root.key = 'Root';
        root.type = 'Boolean';
        root.default = 1;
    option.root = root;
    
        thr.key = 'Threshold';
        thr.type = 'Integer';
        thr.default = NaN;
    option.thr = thr;

        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [.05 .5];
    option.frame = frame;

specif.option = option;
specif.nochunk = 1;
%specif.eachchunk = @eachchunk;
%specif.combineframes = @combineframes;
%specif.afterchunk = @afterchunk;

varargout = mirfunction(@mirlowenergy,x,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
if option.asr
    option.root = 0;
end
if isamir(x,'miraudio')
    if isframed(x)
        x = mirrms(x,'Root',option.root);
    else
        x = mirrms(x,'Frame',option.frame.length.val,option.frame.length.unit,...
                             option.frame.hop.val,option.frame.hop.unit,...
                             'Root',option.root);
    end
end
type = 'mirscalar';


function e = main(r,option,postoption)
if iscell(r)
    r = r{1};
end
if isnan(option.thr)
    if option.asr
        option.thr = .5;
    else
        option.thr = 1;
    end
end
v = mircompute(@algo,get(r,'Data'),option.thr);
fp = mircompute(@noframe,get(r,'FramePos'));
e = mirscalar(r,'Data',v,'Title','Low energy','Unit','/1','FramePos',fp);
e = {e,r};


function v = algo(d,thr)
v = sum(d < repmat(thr*mean(d,2),[1 size(d,2) 1]));
%if not(option.nbframes)
    v = v / size(d,2);
%end


function fp = noframe(fp)
fp = [fp(1);fp(end)];


%function y = eachchunk(orig,option,missing)
%option.nbframes = 1;
%y = mirlowenergy(orig,option);


%function y = combineframes(old,new)  %% not good: last chunk may have different size
%do = mirgetdata(old);
%dn = mirgetdata(new);
%y = set(old,'ChunkData',do+dn);


%function y = afterchunk(orig,frames,postoption)
%d = get(orig,'Data');
%v = mircompute(@afternorm,d,size(frames,2));
%y = set(orig,'Data',v);


%function e = afternorm(d,length)
%e = d/length;