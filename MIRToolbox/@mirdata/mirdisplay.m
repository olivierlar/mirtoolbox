function mirdisplay(d,varargin)
% MIRDATA/DISPLAY display of a MIR data

disp(' ');
v = d.data;
f = d.sr;
n = d.name;
l = d.label;
p = d.pos;
fp = d.framepos;
pp = d.peak.pos;
pm = d.peak.mode;
ld = length(v);
if isempty(d.attack)
    ap = cell(ld);
else
    ap = d.attack.pos;
end
if isempty(d.release)
    rp = cell(ld);
else
    rp = d.release.pos;
end
if isempty(d.track)
    tp = cell(ld);
    tv = cell(ld);
else
    tp = d.track.pos;
    tv = d.track.val;
end
if ld == 0
    disp('No data.');
else
    for i = 1:length(v)
        if nargin < 2
            va = inputname(1);
        else
            va = varargin{1};
        end
        if isempty(va)
            va = 'ans';
        end
        if length(v)>1
            va = [va,'(',num2str(i),')'];
        end
        if not(isempty(l)) && iscell(l) && not(isempty(l{i}))
            lab = ' with label ';
            if isnumeric(l{i})
                lab = [lab,num2str(l{i})];
            else
                lab = [lab,l{i}];
            end
        else
            lab = '';
        end
        disp([va,' is the ',d.title,' related to ',n{i},lab,...
            ', of sampling rate ',num2str(f{i}),' Hz.'])
        if size(v{i},2) == 0
            if isempty(d.init)
                disp('It does not contain any data.');
            else
                disp('It has not been loaded yet.');
            end
        else
            if iscell(d.channels)
                cha = d.channels{i};
            else
                cha = [];
            end
            flag = displot(p{i},v{i},d.abs,d.ord,d.title,fp{i},pp{i},tp{i},tv{i},...
                cha,d.multidata,pm{i},ap{i},rp{i},d.clusters{i});
            if flag
                fig = get(0,'CurrentFigure');
                disp(['Its content is displayed in Figure ',num2str(fig),'.']);
            else
                disp('It does not contain any data.');
            end
        end
    end
end
disp(' ');