function varargout = miremotion(orig,varargin)
% Predicts emotion, simply...

specif.defaultframelength = 1;
specif.defaultframehop = .5;

varargout = mirfunction(@miremotion,orig,varargin,nargout,specif,@init,@main);


function [x type] = init(x,option)
activity = mirrms(x);
valence = mirmode(x);
x = {activity valence};
type = {'miremotion','mirscalar','mirscalar'};


function e = main(x,option,postoption)
rms = get(x{1},'Data');
mode = get(x{2},'Data');
emotion.activity = mircompute(@activity_model,rms);
emotion.valence = mircompute(@valence_model,mode);
e = class(emotion,'miremotion',mirdata(x{1}));
e = set(e,'Title','Emotion');


function v = valence_model(m)
v = m/2+.5;


function a = activity_model(rms)
a = rms*2;