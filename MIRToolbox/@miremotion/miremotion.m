function varargout = miremotion(orig,varargin)
% Predicts emotion along three dimensions and five basic concepts.
% Optional parameters:
%   miremotion(...,'Dimensions',0) excludes all three dimensions.
%   miremotion(...,'Dimensions',3) includes all three dimensions (default).
%   miremotion(...,'Activity') includes the 'Activity' dimension. 
%   miremotion(...,'Valence') includes the 'Valence' dimension. 
%   miremotion(...,'Tension') includes the 'Tension' dimension. 
%   miremotion(...,'Dimensions',2) includes 'Activity' and 'Valence'.
%   miremotion(...,'Arousal') includes 'Activity' and 'Tension'.
%   miremotion(...,'Concepts',0) excludes all five concepts.
%   miremotion(...,'Concepts') includes all five concepts (default).
%   miremotion(...,'Happy') includes the 'Happy' concept.
%   miremotion(...,'Sad') includes the 'Sad' concept.
%   miremotion(...,'Tender') includes the 'Tender' concept.
%   miremotion(...,'Anger') includes the 'Anger' concept.
%   miremotion(...,'Fear') includes the 'Fear' concept.
%   miremotion(...,'Frame',...) predict emotion frame by frame.
%
% Selection of features and coefficients are taken from a study: 
%         Eerola, T., Lartillot, O., and Toiviainen, P. 
%            (2009). Prediction of multidimensional emotional ratings in 
%            music from audio using multivariate regression models. 
%            In Proceedings of 10th International Conference on Music Information Retrieval 
%            (ISMIR 2009), pages 621-626.
%
% The implemented models are based on multiple linear regression with 5 best
% predictors (MLR option in the paper). The box-cox transformations have now been 
% removed until the normalization values have been established with a large sample of music.
% 
% TODO: Revision of coefficients to (a) force the output range between 0 - 1 and 
%    (b) to be based on alternative models and materials (training sets). 
%
% Updated 03.05.2010 TE
%
        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [2 .5];
    option.frame = frame;

        dim.key = 'Dimensions';
        dim.type = 'Integer';
        dim.default = NaN;
        dim.keydefault = 3;
    option.dim = dim;

        activity.key = 'Activity';
        activity.type = 'Boolean';
        activity.default = NaN;
    option.activity = activity;

        valence.key = 'Valence';
        valence.type = 'Boolean';
        valence.default = NaN;
    option.valence = valence;

        tension.key = 'Tension';
        tension.type = 'Boolean';
        tension.default = NaN;
    option.tension = tension;
    
        arousal.key = 'Arousal';
        arousal.type = 'Boolean';
        arousal.default = NaN;
    option.arousal = arousal;
    
        concepts.key = 'Concepts';
        concepts.type = 'Boolean';
        concepts.default = NaN;
    option.concepts = concepts;

        happy.key = 'Happy';
        happy.type = 'Boolean';
        happy.default = NaN;
    option.happy = happy;

        sad.key = 'Sad';
        sad.type = 'Boolean';
        sad.default = NaN;
    option.sad = sad;

        tender.key = 'Tender';
        tender.type = 'Boolean';
        tender.default = NaN;
    option.tender = tender;

        anger.key = 'Anger';
        anger.type = 'Boolean';
        anger.default = NaN;
    option.anger = anger;

        fear.key = 'Fear';
        fear.type = 'Boolean';
        fear.default = NaN;
    option.fear = fear;
    
specif.option = option;
specif.defaultframelength = 2;
%specif.defaultframehop = .5;

specif.combinechunk = {'Average',@nothing};
specif.extensive = 1;

varargout = mirfunction(@miremotion,orig,varargin,nargout,specif,@init,@main);


%%
function [x type] = init(x,option)

option = process(option);

if option.frame.length.val
    hop = option.frame.hop.val;
    if strcmpi(option.frame.hop.unit,'Hz')
        hop = 1/hop;
        option.frame.hop.unit = 's';
    end
    if strcmpi(option.frame.hop.unit,'s')
        hop = hop*get(x,'Sampling');
    end
    if strcmpi(option.frame.hop.unit,'%')
        hop = hop/100;
        option.frame.hop.unit = '/1';
    end
    if strcmpi(option.frame.hop.unit,'/1')
        hop = hop*option.frame.length.val;
    end
    frames = 0:hop:1000000;
    x = mirsegment(x,frames');
elseif isa(x,'mirdesign')
    x = set(x,'NoChunk',1);
end
rm = mirrms(x,'Frame',.046,.5);

le = 0; %mirlowenergy(rm,'ASR');

o = mironsets(x,'Filterbank',15,'Contrast',0.1);
at = mirattacktime(o);
as = 0; %mirattackslope(o);
ed = 0; %mireventdensity(o,'Option1');

fl = mirfluctuation(x,'Summary');
fp = mirpeaks(fl,'Total',1);
fc = 0; %mircentroid(fl);

tp = 0; %mirtempo(x,'Frame',2,.5,'Autocor','Spectrum');
pc = mirpulseclarity(x,'Frame',2,.5); %%%%%%%%%%% Why 'Frame'?? 

s = mirspectrum(x,'Frame',.046,.5);
sc = mircentroid(s);
ss = mirspread(s);
sr = mirroughness(s);

%ps = mirpitch(x,'Frame',.046,.5,'Tolonen');

c = mirchromagram(x,'Frame','Wrap',0,'Pitch',0);    %%%%%%%%%%%%%%%%%%%% Previous frame size was too small.
cp = mirpeaks(c,'Total',1);
ps = 0;%cp;
ks = mirkeystrength(c);
[k kc] = mirkey(ks);
mo = mirmode(ks);
hc = mirhcdf(c);

se = mirentropy(mirspectrum(x,'Collapsed','Min',40,'Smooth',70,'Frame',1.5,.5)); %%%%%%%%% Why 'Frame'?? 

ns = mirnovelty(mirspectrum(x,'Frame',.1,.5,'Max',5000),'Normal',0);
nt = mirnovelty(mirchromagram(x,'Frame',.2,.25),'Normal',0);    %%%%%%%%%%%%%%%%%%%% Previous frame size was too small.
nr = mirnovelty(mirchromagram(x,'Frame',.2,.25,'Wrap',0),'Normal',0);   %%%%%%%%%%%%%%%%%%%% Previous frame size was too small.



x = {rm,le, at,as,ed, fp,fc, tp,pc, sc,ss,sr, ps, cp,kc,mo,hc, se, ns,nt,nr};

type = {'miremotion','mirscalar','mirscalar',...
                     'mirscalar','mirscalar','mirscalar',...
                     'mirspectrum','mirscalar',...
                     'mirscalar','mirscalar',...
                     'mirscalar','mirscalar','mirscalar',...
                     'mirscalar',...
                     'mirchromagram','mirscalar','mirscalar','mirscalar',...
                     'mirscalar',...
                     'mirscalar','mirscalar','mirscalar'};
                 

%%
function e = main(x,option,postoption)

warning('WARNING IN MIREMOTION: The current model of miremotion is not correctly calibrated with this version of MIRtoolbox (but with version 1.3 only).');

option = process(option);
rm = get(x{1},'Data');
%le = get(x{2},'Data');
at = get(x{3},'Data');
%as = get(x{4},'Data');
%ed = get(x{5},'Data');
%fpp = get(x{6},'PeakPosUnit');
fpv = get(x{6},'PeakVal');
%fc = get(x{7},'Data');
%tp = get(x{8},'Data');
pc = get(x{9},'Data');
sc = get(x{10},'Data');
ss = get(x{11},'Data');
rg = get(x{12},'Data');
%ps = get(x{13},'PeakPosUnit');
cp = get(x{14},'PeakPosUnit');
kc = get(x{15},'Data');
mo = get(x{16},'Data');
hc = get(x{17},'Data');
se = get(x{18},'Data');
ns = get(x{19},'Data');
nt = get(x{20},'Data');
nr = get(x{21},'Data');


e.dim = {};
e.dimdata = mircompute(@initialise,rm);
if option.activity == 1
    [e.dimdata e.activity_fact] = mircompute(@activity,e.dimdata,rm,fpv,sc,ss,se);
    e.dim = [e.dim,'Activity'];
else   
    e.activity_fact = NaN;
end
if option.valence == 1
    [e.dimdata e.valence_fact] = mircompute(@valence,e.dimdata,rm,fpv,kc,mo,ns);
    e.dim = [e.dim,'Valence'];
else
    e.valence_fact = NaN;
end
if option.tension == 1
    [e.dimdata e.tension_fact] = mircompute(@tension,e.dimdata,rm,fpv,kc,hc,nr);
    e.dim = [e.dim,'Tension'];
else
    e.tension_fact = NaN;
end

e.class = {};
e.classdata = mircompute(@initialise,rm);
if option.happy == 1
    [e.classdata e.happy_fact] = mircompute(@happy,e.classdata,fpv,ss,cp,kc,mo);
    e.class = [e.class,'Happy'];
else
    e.happy_fact = NaN;
end
if option.sad == 1
    [e.classdata e.sad_fact] = mircompute(@sad,e.classdata,ss,cp,mo,hc,nt);
    e.class = [e.class,'Sad'];
else
    e.sad_fact = NaN;
end
if option.tender == 1
    [e.classdata e.tender_fact] = mircompute(@tender,e.classdata,sc,rg,kc,hc,ns);
    e.class = [e.class,'Tender'];
else
    e.tender_fact = NaN;
end
if option.anger == 1
    [e.classdata e.anger_fact] = mircompute(@anger,e.classdata,rg,kc,se,nr);
    e.class = [e.class,'Anger'];
else
    e.anger_fact = NaN;
end
if option.fear == 1
    [e.classdata e.fear_fact] = mircompute(@fear,e.classdata,rm,at,fpv,kc,mo);
    e.class = [e.class,'Fear'];
else
    e.fear_fact = NaN;
end

e = class(e,'miremotion',mirdata(x{1}));
e = purgedata(e);
fp = mircompute(@noframe,get(x{1},'FramePos'));
e = set(e,'Title','Emotion','Abs','emotions','Ord','magnitude','FramePos',fp);
      
%%      
function option = process(option)
if option.arousal==1
    option.activity = 1;
    option.tension = 1;
    if isnan(option.dim)
        option.dim = 0; 
    end
end
if option.activity==1 || option.valence==1 || option.tension==1
    if isnan(option.activity)
        option.activity = 0;
    end
    if isnan(option.valence)
        option.valence = 0;
    end
    if isnan(option.tension)
        option.tension = 0;
    end
    if isnan(option.concepts)
        option.concepts = 0;
    end
end
if not(isnan(option.dim)) && option.dim
    if isnan(option.concepts)
        option.concepts = 0;
    end
end
if not(isnan(option.concepts)) && option.concepts
    if isnan(option.dim)
        option.dim = 0;
    end
end
if not(isnan(option.dim))
    switch option.dim 
        case 0
            if isnan(option.activity)
                option.activity = 0;
            end
            if isnan(option.valence)
                option.valence = 0;
            end
            if isnan(option.tension)
                option.tension = 0;
            end
        case 2
            option.activity = 1;
            option.valence = 1;
            if isnan(option.tension)
                option.tension = 0;
            end
        case 3
            option.activity = 1;
            option.valence = 1;
            option.tension = 1;
    end
end
if isnan(option.activity)
    option.activity = 1;
end
if isnan(option.valence)
    option.valence = 1;
end
if isnan(option.tension)
    option.tension = 1;
end
if isnan(option.concepts)
    option.concepts = 1;
end
if option.concepts
    option.happy = 1;
    option.sad = 1;
    option.tender = 1;
    option.anger = 1;
    option.fear = 1;
end
if option.happy==1 || option.sad==1 || option.tender==1 ...
        || option.anger==1 || option.fear==1
    if isnan(option.happy)
        option.happy = 0;
    end
    if isnan(option.sad)
        option.sad = 0;
    end
    if isnan(option.tender)
        option.tender = 0;
    end
    if isnan(option.anger)
        option.anger = 0;
    end
    if isnan(option.fear)
        option.fear = 0;
    end
end


%%
function e = initialise(rm)
e = [];

      
function [e af] = activity(e,rm,fpv,sc,ss,se) % without the box-cox transformation, revised coefficients
af = zeros(5,1);

% In the code below, removal of nan values added by Ming-Hsu Chang
af(1) = 0.6664* ((mean(rm(~isnan(rm))) - 0.0559)/0.0337);
tmp = fpv{1};
af(2) =  0.6099 * ((mean(tmp(~isnan(tmp))) - 13270.1836)/10790.655);
tmp = cell2mat(sc);
af(3) = 0.4486*((mean(tmp(~isnan(tmp))) - 1677.7)./570.34);
tmp = cell2mat(ss);
af(4) = -0.4639*((mean(tmp(~isnan(tmp))) - (250.5574*22.88))./(205.3147*22.88)); % New normalisation proposed by Ming-Hsu Chang
af(5) = 0.7056*((mean(se(~isnan(se))) - 0.954)./0.0258);

af(isnan(af)) = [];
e(end+1,:) = sum(af)+5.4861;

function [e vf] = valence(e,rm,fpv,kc,mo,ns) % without the box-cox transformation, revised coefficients
vf = zeros(5,1);
vf(1) = -0.3161 * ((std(rm) - 0.024254)./0.015667);
vf(2) =  0.6099 * ((mean(fpv{1}) - 13270.1836)/10790.655);
vf(3) = 0.8802 * ((mean(kc) - 0.5123)./0.091953);
vf(4) = 0.4565 * ((mean(mo) - -0.0019958)./0.048664);
ns(isnan(ns)) = [];
vf(5) = 0.4015 * ((mean(ns) - 131.9503)./47.6463);
vf(isnan(vf)) = [];
e(end+1,:) = sum(vf)+5.2749;

function [e tf] = tension(e,rm,fpv,kc,hc,nr)
tf = zeros(5,1);
tf(1) = 0.5382 * ((std(rm) - 0.024254)./0.015667);
tf(2) =  -0.5406 * ((mean(fpv{1}) - 13270.1836)/10790.655);
tf(3) = -0.6808 * ((mean(kc) - 0.5124)./0.092);
tf(4) = 0.8629 * ((mean(hc) - 0.2962)./0.0459);
tf(5) = -0.5958 * ((mean(nr) - 71.8426)./46.9246);
tf(isnan(tf)) = [];
e(end+1,:) = sum(tf)+5.4679;


% BASIC EMOTION PREDICTORS

function [e ha_f] = happy(e,fpv,ss,cp,kc,mo)
ha_f = zeros(5,1);
ha_f(1) = 0.7438*((mean(cell2mat(fpv)) - 13270.1836)./10790.655);
ha_f(2) = -0.3965*((mean(cell2mat(ss)) - 250.5574)./205.3147);
ha_f(3) = 0.4047*((std(cell2mat(cp)) - 8.5321)./2.5899);
ha_f(4) = 0.7780*((mean(kc) - 0.5124)./0.092);
ha_f(5) = 0.6220*((mean(mo) - -0.002)./0.0487);
ha_f(isnan(ha_f)) = [];
e(end+1,:) = sum(ha_f)+2.6166;

function [e sa_f] = sad(e,ss,cp,mo,hc,nt)
sa_f = zeros(5,1);
sa_f(1) = 0.4324*((mean(cell2mat(ss)) - 250.5574)./205.3147);
sa_f(2) = -0.3137*((std(cell2mat(cp)) - 8.5321)./2.5899);
sa_f(3) = -0.5201*((mean(mo) - -0.0020)./0.0487);
sa_f(4) = -0.6017*((mean(hc) - 0.2962)./0.0459);
sa_f(5) = 0.4493*((mean(nt) - 42.2022)./36.7782);
sa_f(isnan(sa_f)) = [];
e(end+1,:) = sum(sa_f)+2.9756;

function [e te_f] = tender(e,sc,rg,kc,hc,ns)
te_f = zeros(5,1);
te_f(1) = -0.2709*((mean(cell2mat(sc)) - 1677.7106)./570.3432);
te_f(2) = -0.4904*((std(rg) - 85.9387)./106.0767);
te_f(3) = 0.5192*((mean(kc) - 0.5124)./0.0920);
te_f(4) = -0.3995*((mean(hc) - 0.2962)./0.0459);
te_f(5) = 0.3391*((mean(ns) - 131.9503)./47.6463);
te_f(isnan(te_f)) = [];
e(end+1,:) = sum(te_f)+2.9756;

function [e an_f] = anger(e,rg,kc,se,nr) % 
an_f = zeros(5,1);
%an_f(1) = -0.2353*((mean(pc) - 0.1462)./.1113);
an_f(2) = 0.5517*((mean(rg) - 85.9387)./106.0767);
an_f(3) = -.5802*((mean(kc) - 0.5124)./0.092);
an_f(4) = .2821*((mean(se) - 0.954)./0.0258);
an_f(5) = -.2971*((mean(nr) - 71.8426)./46.9246);
an_f(isnan(an_f)) = [];
e(end+1,:) = sum(an_f)+1.9767;

function [e fe_f] = fear(e,rm,at,fpv,kc,mo)
fe_f = zeros(5,1);
fe_f(1) = 0.4069*((std(rm) - 0.0243)./0.0157);
fe_f(2) = -0.6388*((mean(at) - 0.0707)./0.015689218536423);
fe_f(3) = -0.2538*((mean(cell2mat(fpv)) - 13270.1836)./10790.655);
fe_f(4) = -0.9860*((mean(kc) - 0.5124)./0.0920);
fe_f(5) = -0.3144*((mean(mo) - -0.0019958)./0.048663550639094);
fe_f(isnan(fe_f)) = [];
e(end+1,:) = sum(fe_f)+2.7847;

function fp = noframe(fp)
fp = [fp(1);fp(end)];