function varargout = miremotion(orig,varargin)
% Predicts emotion.

        frame.key = 'Frame';
        frame.type = 'Integer';
        frame.number = 2;
        frame.default = [1 .5];
    option.frame = frame;

specif.option = option;
specif.defaultframelength = 1;
specif.defaultframehop = .5;

varargout = mirfunction(@miremotion,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
x = mirsegment(x,0:5:1000);
[unused,kc] = mirkey(x,'Frame');
md = mirmode(x,'Frame');
ed = mireventdensity(x,'Frame');
sf = mirspectrum(x,'Frame',.1,.5,'Max',5000);
sn = mirnovelty(sf,'KernelSize',16,'Normal',0);
tp = mirtempo(x,'Frame');
br = mirbrightness(x,'Frame','CutOff',110);
rg = mirroughness(x,'Frame');
sf = mirflux(x,'Frame');  %<<<< compare with previous sf
se = mirentropy(x,'Frame');
ose = mirentropy(mirspectrum(x,'Collapsed','Min',40,'Smooth',70));
ar = mirlowenergy(x,'ASR');
fpm = mirpeaks(mirfluctuation(x,'Summary'),'Total',1);
pc = mirpulseclarity(x,'Frame',2,.5);
x = {kc md ed sn tp br rg sf se ose ar fpm pc};
type = {'miremotion','mirscalar','mirscalar','mirscalar','mirscalar',...
        'mirscalar','mirscalar','mirscalar','mirscalar','mirscalar',...
        'mirscalar','mirscalar','mirscalar','mirscalar'};


function e = main(x,option,postoption)
kc = get(x{1},'Data');
md = get(x{2},'Data');
ed = get(x{3},'Data');
sn = get(x{4},'Data');
tp = get(x{5},'Data');
br = get(x{6},'Data');
rg = get(x{7},'Data');
sf = get(x{8},'Data');
se = get(x{9},'Data');
ose = get(x{10},'Data');
ar = get(x{11},'Data');
fpm = get(x{12},'PeakVal');
pc = get(x{13},'Data');
[e.valence e.valence_fact] = ...
                        mircompute(@valence_model,kc,md,ed,sn,tp,br,rg,0);
[e.activity e.activity_fact] = ...
                        mircompute(@activity_model,sf,se,ose,ar,sn,fpm,pc,0);
e = class(e,'miremotion',mirdata(x{1}));
e = set(e,'Title','Emotion');


function [v vf] = valence_model(kc,md,ed,sn,tp,br,rg,ih)
vf(1) = .7135* mean(kc);
vf(2) = .6294*(mean(md)/.382/2+.5);
vf(3) = .5247*mean(ed)/15;
sn(isnan(sn)) = 0;
vf(4) = -.5479*(mean(sn)+440)/855;
vf(5) = .2712*(mean(cell2mat(tp))-40)/160;
vf(6) = -.2703*std(br);
vf(7) = -.2699*mean(rg)/40000;
vf(8) = -.2448*std(ih);
v = sum(vf);


function [a af] = activity_model(sf,se,ose,ar,sn,fpm,pc,ed)
af(1) = .5315*mean(sf)/1000;
af(2) = .6232*mean(se);
af(3) = .4477*ose;
af(4) =  .3833*ar/.71;
sn(isnan(sn)) = 0;
af(5) = .2871*(mean(sn)+440)/855;
af(6) = .2260*(mean(cell2mat(fpm))-3650)/45000;
af(7) = .1877*mean(pc)/.59;
af(8) = .1947*ed;
a = sum(af);