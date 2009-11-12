function varargout = miremotion(orig,varargin)
% Predicts emotion.

        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [0 0];
        frame.keydefault = [1 0];
    option.frame = frame;

specif.option = option;
specif.defaultframelength = 1;
%specif.defaultframehop = .5;

specif.combinechunk = {'Average',@nothing};
specif.framedchunk = 0;

varargout = mirfunction(@miremotion,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)

if option.frame.length.val
    x = mirsegment(x,0:option.frame.length.val:1000000);
end

rm = mirrms(x,'Frame',.046,.5);
le = 0; %mirlowenergy(rm,'ASR');

%o = mironsets(x,'Filterbank',15,'Contrast',0.1);
at = 0; %mirattacktime(o);
as = 0; %mirattackslope(o);
ed = 0; %mireventdensity(o,'Option1');

fl = mirfluctuation(x,'Summary');
fp = mirpeaks(fl,'Total',1);
fc = 0; %mircentroid(fl);

tp = 0; %mirtempo(x,'Frame',2,.5,'Autocor','Spectrum');
pc = 0; %mirpulseclarity(x,'Frame',2,.5);

s = mirspectrum(x,'Frame',.046,.5);
sc = mircentroid(s);
ss = 0; %mirspread(s);
sr = mirroughness(s);

ps = 0; %mirpitch(x,'Frame',.046,.5,'Tolonen');

c = mirchromagram(x,'Frame',.046,.5,'Wrap',0,'Pitch',0);
cp = mirpeaks(c,'Total',1);
ks = mirkeystrength(c);
[k kc] = mirkey(ks);
mo = mirmode(ks);
hc = mirhcdf(c);

se = 0; %mirentropy(mirspectrum(x,'Collapsed','Min',40,'Smooth',70,'Frame',1.5,.5));

ns = 0; %mirnovelty(mirspectrum(x,'Frame',.1,.5,'Max',5000),'Normal',0);
nt = mirnovelty(mirchromagram(x,'Frame',.1,.5),'Normal',0);
nr = 0; %mirnovelty(mirchromagram(x,'Frame',.1,.5,'Wrap',0),'Normal',0);


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
                 

function e = main(x,option,postoption)
rm = get(x{1},'Data');
%le = get(x{2},'Data');
%at = get(x{3},'Data');
%as = get(x{4},'Data');
%ed = get(x{5},'Data');
%fpp = get(x{6},'PeakPosUnit');
fpv = get(x{6},'PeakVal');
%fc = get(x{7},'Data');
%tp = get(x{8},'Data');
%pc = get(x{9},'Data');
sc = get(x{10},'Data');
%ss = get(x{11},'Data');
rg = get(x{12},'Data');
%ps = get(x{13},'Data');
cp = get(x{14},'PeakPosUnit');
kc = get(x{15},'Data');
mo = get(x{16},'Data');
hc = get(x{17},'Data');
%se = get(x{18},'Data');
%ns = get(x{19},'Data');
nt = get(x{20},'Data');
%nr = get(x{21},'Data');

[de e.activity_fact e.valence_fact e.tension_fact] = ...
      mircompute(@models,rm,fpv,sc,rg,cp,kc,mo,hc,nt);
                   
e = class(e,'miremotion',mirdata(x{1}));
e = purgedata(e);
fp = mircompute(@noframe,get(x{1},'FramePos'));
e = set(e,'Title','Emotion','Abs','Concepts','Ord','magnitude',...
          'Data',de,'Concepts',{{{'Activity','Valence','Tension'}'}},...
          'FramePos',fp);

function e = models(rm,fpv,sc,rg,cp,kc,mo,hc,nt)
af(1) = 0.40099*((boxcox(mean(fpv{1}),-.8,1.1011e+04)-2.3546e+07)/5.4559e+03); % NORMALIZATION NEEDS TO BE CALCULATED FROM BOX-COX SCORES
af(2) = 0.25906*((mean(cell2mat(sc)) - 1677.7)./570.34);
af(3) = 1.14410*((mean(rg) - 197.39)./112.72);
af(4) = 0.35787*((mean(cell2mat(cp)) - 8.5321)./2.5899);
af(5) = 0.61027*((mean(hc) - 0.29615)./0.045898);
af(isnan(af)) = [];
e(1,:) = sum(af)+5.4861;

vf(1) = -0.31258 * ((std(rm) - 0.024254)./0.015667);
vf(2) =  0.64626 * ((boxcox(mean(fpv{1}),-.8,1.1011e+04)-2.3546e+07)/5.4559e+03);
vf(3) = 0.88766 * ((mean(kc) - 0.5123)./0.091953);
vf(4) = 0.37939 * ((mean(mo) - -0.0019958)./0.048664);
nt(isnan(nt)) = [];
vf(5) = 0.42506 * ((mean(nt) - 172.84)./41.334);
vf(isnan(vf)) = [];
e(2,:) = sum(vf)+5.2749;

tf(1) = 0.6333 * ((std(rm) - 0.024254)./0.015667);
tf(2) =  -0.42467 * ((boxcox(mean(fpv{1}),-.8,1.1011e+04)-2.3546e+07)/5.4559e+03);
tf(3) = 0.55581 * ((mean(cell2mat(sc)) - 1677.7)./570.34);
tf(4) = -0.82602 * ((mean(kc) - 0.5123)./0.091953);
tf(5) = -0.87464 * ((mean(nt) - 172.84)./41.334);
tf(isnan(tf)) = [];
e(3,:) = sum(tf)+5.4679;

e = {e af vf tf};


function v=boxcox(x,lamda,xdot)
%Syntax: v=boxcox(x,lamda,xdot)
%______________________________
%
% Makes the Box-Cox transformation of a data set x.
%
% v is the transformed data vector.
% x is the data set.
% lamda is the parameter of the transformation.
% xdot is the geometric mean of the data.
%
% Alexandros Leontitsis
% Institute of Mathematics and Statistics
% University of Kent at Canterbury
% Canterbury
% Kent, CT2 7NF
% U.K.
% University e-mail: al10@ukc.ac.uk (until December 2001)
% Lifetime e-mail: leoaleq@yahoo.com
% Homepage: http://www.geocities.com/CapeCanaveral/Lab/1421
%
% June 15, 2001.

if nargin<1 | isempty(x)==1
   error('You shoud provide a data set.');
else
   % x must be a vector
   if min(size(x))>1
      error('Invalid data set.');
   end
   x=x(:);
   % n is the data set length
   n=length(x);
end

if nargin<2 | isempty(lamda)==1
   lamda=0;
else
   %lamda must be either a scalar or a vector
   if min(size(lamda))>1
      error('lamda must be either a scalar or a vector.');
   end
end

if nargin<3 | isempty(xdot)==1
   xdot=1;
else
   % xdot must be a scalar
   if sum(size(xdot))>2
      error('xdot must be a scalar.');
   end
end

% Box-Cox transformation
for i=1:length(lamda)
   if lamda(i)~=0
      v(:,i)=(x.^lamda(i)-1)/(lamda(i)*xdot^(lamda(i)-1));
   else
      v(:,i)=xdot*log(x);
   end
end


function fp = noframe(fp)
fp = [fp(1);fp(end)];