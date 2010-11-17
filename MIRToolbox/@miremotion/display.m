function display(e)
% MIREMOTION/DISPLAY display of emotion contents
d = get(e,'Dim');
dd = get(e,'DimData');
c = get(e,'Class');
cd = get(e,'ClassData');
fp = get(e,'FramePos');
n = get(e,'Name');
if not(isempty(d))
    s = mirscalar(e);
    s = set(s,'Pos',repmat({{d'}},size(dd)),'Data',dd,...
        'Title','Dimensional Emotion','Abs','Emotion Dimensions')
end
if not(isempty(c))
    s = mirscalar(e);
    s = set(s,'Pos',repmat({{c'}},size(cd)),'Data',cd,...
        'Title','Basic Emotion Set','Abs','Emotion Classes')
end