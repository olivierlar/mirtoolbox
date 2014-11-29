function display(s,ax,songs,suffix)
% SCALAR/DISPLAY display the values of a scalar object

ST = dbstack;
if strcmp(ST(end).file,'arrayviewfunc.m')
    disp('To display its content in a figure, evaluate this variable directly in the Command Window.');
    return
end

disp(' ');
p = get(s,'Pos');
v = get(s,'Data');
m = get(s,'Mode');
t = get(s,'Title');
n = get(s,'Name');
n2 = get(s,'Name2');
u = get(s,'Unit');
fp = get(s,'FramePos');
pp = get(s,'PeakPos');
leg = get(s,'Legend');
legm = get(s,'MultiData');
cha = get(s,'Channels');

if nargin<3 || isempty(songs)
    songs=1:length(v);
end

for song = 1:length(songs)  %For each audio file
    i=songs(song);

    vi = v{i};
    if isempty(m)
        mi = [];
    else
        mi = m{i};
    end
    fpi = fp{i};
    if isempty(pp)
        ppi = [];
    else
        ppi = pp{i};
    end
    
    if isempty(vi)
        disp(['The ',t,' related to file ',n{i},...
            ' does not contain any data.']);
        return
    end
        
    if size(vi,1) == 1 && size(vi,2) == 1 && size(vi,3) == 1 && ...
            (not(iscell(vi)) || (size(vi{1},1) == 1 && ...
            size(vi{1},2) == 1 && size(vi{1},3) == 1)) && ...
            (not(iscell(vi{1})) || (size(vi{1}{1},1) == 1 && ...
            size(vi{1}{1},2) == 1 && size(vi{1},3) == 1))
        % Simple results, returned directly in the Command Window
        
        if iscell(vi)
            vi = vi{1}; % There is only one segment, so let's look at it
            if not(isempty(mi))
                mi = mi{1};
            end
        end
        if iscell(vi)
            if size(vi,2) > 1
                error('BUG IN DISPLAY'); %test if this condition exists
            end
            vi = vi{1}; % There is only one frame, so let's look at it
            if not(isempty(mi))
                mi = mi{1};
            end
        end
        if size(vi,1) == 1 && size(vi,3) == 1
            if isempty(leg)
                r = num2str(vi);
            else
                r = leg{vi};
            end
            if not(isempty(mi))
                u = legm{mi};
            end
            if strcmp(u,'/1')
                u = '';
            end
            if strcmp(r,'NaN')
                r = 'undefined';
                u = '';
            end
            if not(isempty(n2))
                disp(['The ',t,' between files ',n{1},' and ',n2{i},...
                    ' is ',r,' ',u]);
            else
                disp(['The ',t,' related to file ',n{i},' is ',r,' ',u]);
            end
        else
            vi = squeeze(vi);
            disp(['The ',t,' related to file ',n{i},' are:']);
            for j = 1:size(vi,1)
                if isempty(leg)
                    r = num2str(vi(j));
                else
                    r = leg{vi(j)};
                end
                if not(isempty(mi))
                    u = legm{mi(j)};
                end
                if strcmp(r,'NaN')
                    r = 'undefined';
                    u = '';
                end
                disp(['Value #',num2str(j),': ',r,' ',u]);
            end
        end
        
    else
        %Graphical display
        
        if nargin<2 || isempty(ax) || ischar(ax)
            figure
        else
            axes(ax)
        end
        
        if not(iscell(vi))
            vi = {vi};
            fpi = {fpi};
            ppi = {ppi};
        end
        
        nl = 0;
        for j = 1:length(vi)
            if 0 %iscell(vi{j}) && length(vi{j}) == 1
                vi{j} = vi{j}{1};
            end
            if ~isempty(vi{j})
                nl = size(vi{j},1);    % Number of bins
                nc = size(vi{j},2);    % Number of frames
                l = size(vi{j},3);     % Number of channels
                il = (1-0.15)/l;
                break
            end
        end
        if ~nl
            break
        end
        
        if nc==1 && l>1        % If one frame and several channels
            xlab = 'Channels'; % channels will be represented in the x axis
            l = 1;
        else
            xlab = 'Temporal location of events (in s.)';
        end
        
        varpeaks = 0; % Variable data size over frames?
        if ~isempty(mi)
            varpeaks = 1;
        elseif iscell(vi{1})
            for j = 1:length(vi)
                for k = 1:l
                    for h = 1:size(vi{j},2)
                        if length(vi{j}{1,h,k}) > 1
                            varpeaks = 1;
                        end
                    end
                end
            end
            if not(varpeaks)
                for j = 1:length(vi)
                    vj = zeros(size(vi{j}));
                    for k = 1:l
                        for h = 1:size(vi{j},2)
                            if isempty(vi{j}{1,h,k})
                                vj(1,h,k) = NaN;
                            else
                                vj(1,h,k) = vi{j}{1,h,k};
                            end
                        end
                    end
                    vi{j} = vj;
                end
            end
        end
        
        if varpeaks         % Peaks displayed with diamonds
            diamond = 1;
            set(gca,'NextPlot','replacechildren')
            hold on
        else
            diamond = 0;
            hold all
        end
        
        for k = 1:l         % For each channel
            if l>1
                subplot('Position',[0.1 (k-1)*il+0.1 0.89 il-0.02])
            end
            hold on
            
            vold = NaN;     % Buffers used to link curves between segments
            fold = NaN;
            
            for j = 1:length(vi)     % for each segment
                vj = vi{j};
                
                fpj = fpi{j};
                if isempty(fpj)
                    continue
                end
                if isempty(ppi) || length(ppi)<j
                    ppj = [];
                else
                    ppj = ppi{j};
                end
                
                if strcmp(xlab,'Channels')
                    %vj = squeeze(vj); % does not work with lowenergy
                    %ppj = squeeze(ppj);
                    vj = shiftdim(vj,1);
                    ppj = shiftdim(ppj,1);
                    fpj = cha{i};
                end
                
                minvi = Inf;    % The limits of the rectangle bordering segments
                maxvi = -Inf;
                
                if diamond
                    % Variable data size for each frame
                    legobj = [];
                    legval = [];
                    
                    if not(strcmp(xlab,'Channels'))
                        if iscell(fpi(1,j))
                            fpj = fpi{j};
                        else
                            error('BUG IN DISPLAY'); %test if this condition exists
                            fpj = fpi(:,j);
                        end
                    end
                    
                    for h = 1:size(vj,2)    % for each frame
                        if not(isempty(vj{1,h,k}))
                            if isempty(mi)
                                % No use of legends
                                if not(isempty(vj{1,h,k}))
                                    plot(mean(fpj(:,h),1),vj{1,h,k},'d',...
                                        'LineWidth',2,'MarkerFaceColor','b');
                                end
                                
                            else
                                % Use of legends
                                mj = mi{j}{1,h,k};
                                lmj = length(mj);
                                for kk = 1:lmj
                                    handle = plot(mean(fpj(:,h),1),vj{1,h,k}(kk),'d',...
                                        'Color',num2col(mj(kk)),...
                                        'MarkerSize',7+(lmj-kk)*2,...
                                        'MarkerFaceColor',num2col(mj(kk)));
                                    if isempty(find(legval == mj(kk)))
                                        legval = [legval mj(kk)];
                                        legobj = [legobj handle];
                                    end
                                end
                            end
                            
                            %if not(isempty(ppi)) ...    % Peaks display
                            %&& not(isempty(find(ppi{1} == j)))
                            %    plot(mean(fpj(:,h),1),vj{1,h,k},'+r')
                            %end
                        end
                        
                        if not(isempty(vj{1,h,k}))
                            minvi = min(minvi,min(vj{1,h,k}));
                            maxvi = max(maxvi,max(vj{1,h,k}));
                        end
                    end
                    
                    if (exist('legm') == 1) && not(isempty(legm))
                        legend(legobj,legm{legval},'Location','Best')
                    end
                    
                else
                    % Constant data size along frames
                    
                    for h = 1:nl    % For each dimension represented
                        %if not(isempty(vj(h,:,k)))
                        if isnan(vold)
                            plot(mean(fpj,1),vj(h,:,k)',...
                                '-','Color',num2col(h))
                        else
                            plot([fold mean(fpj,1)],[vold(h) vj(h,:,k)]',...
                                '+-','Color',num2col(h))
                            % Current curve linked with curve from
                            % previous segment
                        end
                        if h == nl
                            if isempty(vj)
                                vold = NaN;
                                fold = NaN;
                            else
                                vold = vj(:,end,k);
                                fold = mean(fpj(:,end),1);
                            end
                        end
                        % else
                        %     vold = NaN;
                        %     fold = NaN;
                        % end
                    end
                    
                    if isa(s,'mirpitch')
                        ps = get(s,'Start');
                        pe = get(s,'End');
                        pm = get(s,'Mean');
                        deg = get(s,'Degrees');
                        if length(ps)>=i && ~isempty(ps{i})              % Note display
                            for h = 1:nl;
                                psh = ps{i}{j}{1,h,k};
                                peh = pe{i}{j}{1,h,k};
                                pmh = pm{i}{j}{1,h,k};
                                plot([fpj(1,psh);fpj(2,peh)],...
                                     [pmh;pmh],'rd-','LineWidth',2)
                                if ~isempty(deg) && ~isempty(deg{i})
                                    for hh = 1:length(pmh)
                                        text(mean(fpj(:,psh(hh))),pmh(hh)+40,...
                                             num2str(deg{i}{j}{1,h,k}(hh)),...
                                             'FontSize',15,'Color','r');
                                    end
                                end
                            end
                        end
                    end
                    if not(isempty(ppj))                % Peaks display
                        for h = 1:nl;
                            pph = ppj{1,h,k};
                            plot(mean(fpj(:,pph)),vj(h,pph,k),'or')
                        end
                    end
                    
                    minvi = min(minvi,min(min(vj(h,:,k))));
                    maxvi = max(maxvi,max(max(vj(h,:,k))));
                end
                
                if length(vi)>1 && size(vj,2)>1 && minvi<Inf && maxvi>-Inf
                    % Display of the segment rectangle
                    rectangle('Position',[fpj(1),minvi,...
                        fpj(end)-fpj(1),maxvi-minvi]+1e-16,...
                        'EdgeColor',num2col(j),'Curvature',.1,'LineWidth',1)
                end
            end
            if k == l
                title([t,', ',n{i}])
            end
            if k == 1
                xlabel(xlab)
            end
            if not(isempty(leg))
                set(gca,'ytick',(1:length(leg)))
                set(gca,'yticklabel',leg);
            end
            if l > 1
                pos = get(gca,'Position');
                hfig = axes('Position',[pos(1)-.05 pos(2)+pos(4)/2 .01 .01],...
                    'Visible','off');
                text(0,0,num2str(cha{i}(k)),'FontSize',12,'Color','r')
            end
        end
        if i == 1
            if strcmp(u,'/1')
                u = ' (between 0 and 1)';
            elseif isempty(u)
                u = '';
            else
                u = [' (in ',u,')'];
            end
        end
        ylabel(['coefficient value', u])
        if nl>1
            legnd = cell(nl,1);
            for j = 1:nl
                if isempty(p) || isempty(p{1})
                    legnd{j} = num2str(j);
                else
                    legnd{j} = p{i}{1}{j};
                end
            end
            legend(legnd,'Location','Best')
        end
        fig = get(0,'CurrentFigure');
        if isa(fig,'matlab.ui.Figure')
            fig = fig.Number;
        end
        disp(['The ',t,' related to file ',n{i},...
            ' is displayed in Figure ',num2str(fig),'.']);
        if nargin>3
            saveas(fig,[n{i},suffix],'psc2');
            disp(['and is saved in file ',n{i},suffix]);
        end
        axis tight
        axis 'auto y'
    end
    if nargin>1 && ischar(ax)
        print(ax)
    end
end
disp(' ');
drawnow