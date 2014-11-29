function mirdisplay(d,va,axis,songs,suffix)
% MIRDATA/DISPLAY display of a MIR data

disp(' ');
v = d.data;
f = d.sr;
n = d.name;
l = d.label;
p = d.pos;
fp = d.framepos;
ld = length(v);
pp = d.peak.pos;
pm = d.peak.mode;
if isempty(pp)
    pp = cell(ld);
    pm = cell(ld);
end
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
    if nargin<4 || isempty(songs)
        songs=1:length(v);
    end
    
    for song = 1:length(songs)  %For each audio file
        i=songs(song);
        
        if nargin < 2
            va = inputname(1);
        end
        if isempty(va)
            va = 'ans';
        end
        if length(v)>1
            vai = [va,'(',num2str(i),')'];
        else
            vai = va;
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
        disp([vai,' is the ',d.title,' related to ',n{i},lab,...
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
                cha,d.multidata,pm{i},ap{i},rp{i},d.clusters{i},axis);
            if flag
                fig = get(0,'CurrentFigure');
                if isa(fig,'matlab.ui.Figure')
                    fig = fig.Number;
                end
                disp(['Its content is displayed in Figure ',num2str(fig),'.']);
                if nargin>4 && ~isempty(suffix)
                    saveas(fig,[n{i},suffix],'psc2');
                    disp(['and is saved in file ',n{i},suffix]);
                end
            else
                disp('It does not contain any data.');
            end
        end
    end
end
disp(' ');