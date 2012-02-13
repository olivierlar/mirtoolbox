function MIRplayer(arg,select)
%Usage:
%MIRplayer(a, varargin)
%where
%a= either session data saved by previous execution of MIRplayer,
%   or struct or mirstruct (features),
%------

%TODO
% implement in OOP
% display seconds of the playback
% display frame length and hop factor of the feature
% possibility to add more axes/delete axes, select an axis on which to put a feature
% save settings (axes, selected features, colors)
% display help in a toolbar
% visualize multidimensional data: mfcc, tempo, chromagram
% print the current figure in a file (good for presentations)



if ischar(arg) %session data given in a file
    if (exist(strcat(arg,'/featureInfo.mat'),'file') ~=2) || ~exist(strcat(arg,'/1.wav'),'file')
        error('Session data file not found in a MATLAB file.');
    end
    disp('Loading the session. Please wait...');
    load(strcat(arg,'/featureInfo.mat'));
    if ~exist('features','var') || ~exist('songNames','var') || ~exist('xlims','var') || ~exist('songSampling','var')
        error('Session feature data file is not compatible with MIRplayer');
    end
    nBins=size(features.distribution,2);
    songs=1:length(songNames);
else
    %the first argument should be the feature set
    if ~isstruct(arg) && ~iscell(arg) && ~isa(arg,'mirscalar')
        error('The first input argument should be struct or mirscalar variable.');
    end
        
    if nargin>1
        songs=int8(unique(select));
        
        %if any(songs)<1 || any(songs)>length(songNames)
        %    error('The second input argument should be an integer array of songs to be included in the analysis.');
        %end
    else
        songs=[];
    end
        
    nBins=100; %resolution of the feature distributions (for visualization)
    smoothingFactor=2;%round(max(3,nBins/10)); %for median filtering the distribution
    features=getFeatureInfo(arg, nBins,smoothingFactor,songs);
    clear('arg');
    
    %TODO: Audio could also be shown without features. Just to be able to check
    %out the songs before extracting features.
    if isempty(features)
        error('Please provide a feature set with at least one mirscalar type of feature.');
    end
    
end


global fig
global Fs
global pointerH
global player
global CurrentSample
global CurrentFrame
global framePos
global xlim
global PauseButton
global PPSstate
global smallPointerH
global sliderH
global aH
global followPointerButton


framediff=0;
song=0;
songImage=0;
pointerY=[];
ylim=[0,1];
CurrentFrameStartPosition=0;

pointOnSlider=[];

mirchunklim=5000000;

nFeatures=length(features.names);

guiColor=[.9,.9,.9];
pointerColor1='k'; %Black
peakColor='r'; %Red
featureColors=[0,0,0;guiColor;0,.6,0;.9,.7,0;0,0,.9;.8,0,0];%0,.5,0;.9,.6,0;0,0,.8;.6,0,0]; %available rgb values for plotting the features
peakColors=[guiColor;0,.5,0;.9,.6,0;0,0,.8;.6,0,0]; %available rgb values for plotting the features
pointerColor=pointerColor1;
pointerAlphaDefault=.2;
peakAlphaDefault=.1;
%downsample a to save memory when plotting audio
downSampleRate=1000; %Could be nice resolution

songColor=[.85,.85,.85];
maxFrameUpdateFrequency=.1; %seconds
zoomFactorDefault=2/3;
%feature=cell(nFeatures,1);
%featureCurvePos=cell(nFeatures,1);
selectFeatureButtons=cell(nFeatures,1);
selectPeaksButtons=cell(nFeatures,1);
mainFeature=0;
selectedFeatures=2*ones(nFeatures,1);
selectedFeatureInd=0;

selectedPeaks=ones(nFeatures,1);

%this is the rugged player implemented in MATLAB
%create play, pause and stop icons for buttons
%s=30;

s=get(0,'ScreenSize'); s=round(s(3)/150)*2; %size of the rectangular area
stopIcon=zeros(s,s,3);
pauseIcon=zeros(s,s,3);
pauseIcon(:,s/2-1:s/2+1,:)=NaN;
tri=triu(ones(s/2,s/2),1);
playIcon=zeros(s,s); playIcon(1:s/2,1:2:s)=tri;
playIcon(1:s/2,2:2:s)=tri;
playIcon=repmat(playIcon+playIcon(s:-1:1,:),[1,1,3]);
playIcon(playIcon==1)=NaN;

%% CREATE GUI

fig    =   figure(...       % the main GUI figure
    'MenuBar','fig', ...
    'Toolbar','none', ...
    'HandleVisibility','callback', ...
    'Units', 'normalized', ...
    'OuterPosition', [0,0,1,1], ...
    'PaperUnits', 'inches', ...
    'PaperOrientation', 'landscape', ...
    'PaperSize', [16,9], ...
    'Name', mfilename, ...
    'NumberTitle','off', ...
    'WindowButtonDownFcn', @startDragFcn, ...
    'Color', guiColor);
menuItems = uimenu('Parent',fig,'Label','File');
%uimenu(menuItems,'Label','Print','Callback','printpreview');
uimenu(menuItems,'Label','Save session','Callback',@saveSession);
uimenu(menuItems,'Label','Quit session','Separator','on','Accelerator','Q','Callback',@quitSession);





    function saveSession(hObject,eventdata)
        folderName=inputdlg('Folder name','Save session');
        if isempty(folderName{1})
            return
        else
            mkdir(folderName{1});
            save(strcat(folderName{1},'/featureInfo'), 'features','songNames','xlims','songSampling');
        end
        
        for i=1:length(songs)
            ind=songs(i);
        %if withMiraudio
        %    wavwrite(songData{ind}{1},songSampling{ind},16,strcat(folderName{1},'/',num2str(i),'.wav'));
        %else
        %    copyfile(strcat(a,'/',num2str(ind),'.wav'), strcat(folderName{1},'/',num2str(i),'.wav'));
        %end
        end
            
        
    end

    function quitSession(hObject,eventdata)
        quitting = questdlg('Save session before quitting?','Quit','Don''t save','Cancel','Save', 'Save');
        if isempty(quitting) || isequal(quitting,'Cancel')
            return;
        elseif isequal(quitting,'Save')
            java.lang.Runtime.getRuntime.gc; %Java garbage collection
            saveSession;
            if ishandle(fig) 
                close(fig); 
            end
        else
            if ishandle(fig)
                close(fig);
            end
        end
    end

mainPanelPos=[.16,.08,.8,.72];
MainPanel = uipanel(...
    'Parent', fig, ...
    'Units', 'normalized', ...
    'Clipping', 'off', ...
    'HandleVisibility', 'callback', ...
    'Position',mainPanelPos, ...
    'BorderType','none', ...
    'BackGroundColor', guiColor, ...
    'Visible','on');
ControlPanel = uipanel(...
    'Parent', fig, ...
    'Units', 'normalized', ...
    'Title', 'SELECT A SONG', ...
    'TitlePosition','centertop', ...
    'FontSize', 10, ...
    'FontUnits', 'normalized', ...
    'Clipping', 'off', ...
    'HandleVisibility', 'callback', ...
    'Position',[.16,.82,.8,.16], ...
    'BorderType','none', ...
    'BackGroundColor', guiColor, ...
    'Visible','on');
FeaturePanel = uipanel(...
    'Parent', fig, ...
    'Title', 'SELECT FEATURES AND PEAKS', ...
    'FontSize', 10, ...
    'FontUnits', 'normalized', ...
    'Units', 'normalized', ...
    'Clipping', 'on', ...
    'HandleVisibility', 'callback', ...
    'Position',[.01, .98-length(features.names)*.035, mainPanelPos(1)-.05, length(features.names)*.035], ...
    'BorderType','none', ...
    'BackGroundColor', guiColor, ...
    'Visible','on');
DistPanel = uipanel(...
    'Parent', fig, ...
    'Title', 'FEATURE DISTRIBUTION', ...
    'FontSize', 10, ...
    'FontUnits', 'normalized', ...
    'Units', 'normalized', ...
    'Clipping', 'off', ...
    'HandleVisibility', 'callback', ...
    'Position',[.01,.08,mainPanelPos(1)-.05,.15], ...
    'BackGroundColor', guiColor, ...
    'BorderType','none', ...
    'Visible','on');
outerPos=get(fig,'OuterPosition');
ratioPos=outerPos(3)*mainPanelPos(3)/(outerPos(4)*mainPanelPos(4));
aH      =   axes(...         % the axes for plotting
    'Parent', MainPanel, ...
    'Units', 'normalized', ...
    'HandleVisibility','callback', ...
    'Position',[0 0 1 1]);
xlabel(aH,'Time (s)');
ylabel(aH,'Feature value');
PPSsize=[.15,.45];

% Create the button group.
PPSbuttons = uibuttongroup('Parent',ControlPanel, ...
    'Position',[0,.4,PPSsize], ...
    'BorderType','none', ...
    'Visible','on', ...
    'BackGroundColor', guiColor, ...
    'SelectionChangeFcn',@selectPPS);

PlayButton  =   uicontrol(...
    'Parent', PPSbuttons, ...
    'Style','Toggle', ...
    'CData',playIcon, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0 0 .3 1], ...%[.1,0.85,buttonSize],...
    'Tag','play');
PauseButton  =   uicontrol(...
    'Parent', PPSbuttons, ...
    'Style','Toggle', ...
    'CData',pauseIcon, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[.35 0 .3 1],...
    'Tag','pause');
StopButton  =   uicontrol(...
    'Parent', PPSbuttons, ...
    'Style','Toggle', ...
    'CData',stopIcon, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[.7 0 .3 1],...
    'Tag','stop');
audioPopupmenuH=   uicontrol(...    % list of available audio
    'Parent', ControlPanel, ...
    'Units','normalized',...
    'Position',[.25 .7 .5 .1],...
    'Callback', @selectSong, ...
    'HandleVisibility','callback', ...
    'String',features.songNames,...
    'TooltipString','Available audio files', ...
    'Style','popupmenu');


%initialize graphics in aH
songH=line( ...
    'Parent',aH, ...
    'XData',[0,0], ...
    'YData',[0,0], ...
    'Color',[1,1,1]);


pointerH=patch( ...
    'Parent',aH, ...
    'XData',[0,0,0,0], ...
    'YData',[0,0,0,0], ...
    'LineWidth',1, ...
    'FaceColor',pointerColor, ...
    'EdgeColor',pointerColor, ...
    'FaceAlpha',0, ...
    'EdgeAlpha',0);

% Create zoom button group.
axesPos=get(aH,'Position');

zoomButtons = uipanel('Parent',ControlPanel, ...
    'Position',[0,0,1,.6], ...
    'BorderType','none', ...
    'BackGroundColor', guiColor, ...
    'Visible','on');

zoomInIcon=load(fullfile(matlabroot, '/toolbox/matlab/icons/zoomplus.mat'));
zoomInIcon=zoomInIcon.cdata;
zoomOutIcon=load(fullfile(matlabroot, '/toolbox/matlab/icons/zoomminus.mat'));
zoomOutIcon=zoomOutIcon.cdata;

zoomInButton  =   uicontrol(...
    'Parent', zoomButtons, ...
    'Style','PushButton', ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[.943 .5 .03 .5], ...%[.1,0.85,buttonSize],...
    'CData',zoomInIcon, ...
    'Tag', 'in', ...
    'CallBack',@zoomAxes);
zoomOutButton  =   uicontrol(...
    'Parent', zoomButtons, ...
    'Style','PushButton', ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[.973 .5 .03 .5], ...%[.1,0.85,buttonSize],...
    'CData',zoomOutIcon, ...
    'Tag', 'out', ...
    'CallBack',@zoomAxes);
%figPos=get(fig,'Position');
sliderAxes  =   axes(...
    'Parent', zoomButtons, ...
    'Units','Normalized',...
    'HandleVisibility','callback', ...
    'XTick',[], ...
    'YTick',[], ...
    'Xlim',[0,1], ...
    'Ylim', [0,1], ...
    'Position',[0,0,1,.5]);%, ...
%'Visible', 'off', ...);
sliderH=patch( ...
    'Parent',sliderAxes, ...
    'XData',[0,1,1,0], ...
    'YData',[0,0,1,1], ...
    'LineWidth',1, ...
    'FaceColor',pointerColor, ...
    'EdgeColor',pointerColor, ...
    'FaceAlpha',.2, ...
    'EdgeAlpha',.2);%, ...
%'Visible', 'off');
followPointerButton  =   uicontrol(...
    'Parent', zoomButtons, ...
    'Style','CheckBox', ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'TooltipString','Follow pointer position', ...
    'Position',[.915 .5 .03 .5]);%, ...
%'CallBack',@followPointer);
songThumbnailH=line( ...
    'Parent',sliderAxes, ...
    'XData',[0,0], ...
    'YData',[0,0], ...
    'Color',[1,1,1]);
smallPointerH=line( ...
    'Parent',sliderAxes, ...
    'XData',[0,0], ...
    'YData',[0,1], ...
    'LineWidth',1, ...
    'Color',[1, .5, .5]);
z=zoom(aH);
setAxesZoomMotion(z,aH,'horizontal');

distAxes=axes(...
    'Parent',DistPanel, ...
    'Units','normalized',...
    'Position', [-.01,0,1,.95], ...
    'Xlim',[0,1], ...
    'Ylim',[0,1], ...
    'Xtick',[], ...
    'Ytick',[], ...
    'LineWidth', .00001, ...
    'Color', guiColor);
featureDistPatch=patch(...
    'Parent',distAxes, ...
    'YData',zeros(1,nBins+2), ...
    'XData',[0:1/(nBins-1):1,1,0], ...
    'FaceAlpha',.2, ...
    'EdgeAlpha',.4, ...
    'FaceColor','r', ...
    'EdgeColor','r');
songDistPatch=patch(...
    'Parent',distAxes, ...
    'YData',[0,0], ...
    'XData',[0:1/(nBins-1):1,1,0], ...
    'FaceAlpha',.2, ...
    'EdgeAlpha',.4, ...
    'FaceColor','g', ...
    'EdgeColor','g');

text('Parent',distAxes, ...
    'String','song', ...
    'FontSize',8, ...
    'FontUnits','normalized',...
    'Units','normalized', ...
    'Position', [0,-.2], ...
    'Color', [0,.5,0]);
text('Parent',distAxes, ...
    'String','all songs', ...
    'FontSize',8, ...
    'FontUnits','normalized',...
    'Units','normalized', ...
    'Position', [.2,-.2], ...
    'Color', [.5,0,0]);
noDataText=text('Parent',distAxes, ...
    'String','NO SONG DATA', ...
    'Units','normalized', ...
    'Position',[.5,.5], ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment','Middle', ...
    'FontSize',12, ...
    'FontUnits', 'normalized', ...
    'FontWeight', 'bold', ...
    'Color',[0,.5,0], ...
    'Visible', 'off');

FeatureActionsH = axes(...
    'Parent', FeaturePanel, ...
    'Units', 'normalized', ...
    'Clipping', 'on', ...
    'HandleVisibility', 'callback', ...
    'Position',[0,0,1,1], ...
    'Xlim',[-.01,10], ...
    'Ylim',[-.01,length(features.names)], ...
    'Xtick',[], ...
    'Ytick',[], ...
    'LineWidth', .00001, ...
    'Color', guiColor, ...
    'Visible','on');

vertRect=[.1,.15;.9,.15;.9,.85;.1,.85];

vrCreate=[0;1;1;0];
colorVertices=[[vertRect(:,1)+1+vrCreate;vertRect(:,1)+3+vrCreate;vertRect(:,1)+5;vertRect(:,1)+6;vertRect(:,1)+7;vertRect(:,1)+8],[repmat(vertRect(:,2),6,1)];1.5,.25;2,.25;2.5,.25;2,.75;3.5,.75;4,.75;4.5,.75;4,.25];
colorFaces=[1,2,3,4;5,6,7,8;9,10,11,12;13,14,15,16;17,18,19,20;21,22,23,24;25,26,27,28;29,30,31,32];
colorCdata=[zeros(4,3);ones(4,3);repmat(featureColors(3,:),4,1);repmat(featureColors(4,:),4,1);repmat(featureColors(5,:),4,1);repmat(featureColors(6,:),4,1);ones(4,3);repmat(guiColor,4,1)];
peakColorVertices=[[vertRect(:,1)+3;vertRect(:,1)+4;vertRect(:,1)+5;vertRect(:,1)+6;vertRect(:,1)+7+vrCreate], [repmat(vertRect(:,2),5,1)];7.5,.25;8,.25;8.5,.25;8,.75];
peakColorFaces=[1,2,3,4;5,6,7,8;9,10,11,12;13,14,15,16;17,18,19,20;21,22,23,24];
peakColorCdata=[repmat(peakColors(2,:),4,1);repmat(peakColors(3,:),4,1);repmat(peakColors(4,:),4,1);repmat(peakColors(5,:),4,1);ones(4,3);repmat(guiColor,4,1)];

colorPatch=patch('Parent',FeatureActionsH, ...
    'Vertices',colorVertices, ...
    'Faces',colorFaces, ...
    'FaceColor','flat','EdgeColor','flat',...
    'FaceVertexCData',colorCdata, ...
    'Visible','off','LineWidth',.01,'ButtonDownFcn', @setFeature);


for i=1:nFeatures
    
    featureH{i}=line( ...
        'Parent',aH, ...
        'XData',[0,0], ...
        'YData',[0,0], ...
        'Color','k');
    
    selectFeatureButtons{i}=patch('Parent',FeatureActionsH,'XData',[.1,.9,.1],'YData',nFeatures-i+[.15, .5,.85],'FaceColor',guiColor,'EdgeColor','k','Tag', int2str(i),'ButtonDownFcn', @selectFeature,'Visible','on');
    
    featureText{i}=text( ...
        'Parent',FeatureActionsH, ...
        'Units','data', ...
        'String',features.names{i}, ...
        'Position',[1.5 nFeatures-i+.15], ...
        'HorizontalAlignment','Left', ...
        'VerticalAlignment','Bottom', ...
        'FontSize',10, ...
        'Color','k', ...
        'Tag', int2str(i), ...
        'ButtonDownFcn', @selectFeature, ...
        'Visible','on');
    
    if features.hasPeaks(i)
        %selectPeaksButtons{i}=patch('Parent',FeatureActionsH,'XData',[9.2,9.8,9.8,9.2],'YData',nFeatures-i+[.2, .2,.8,.8],'FaceColor',guiColor,'EdgeColor','k','EdgeAlpha',.5','Tag', int2str(i),'ButtonDownFcn', @selectPeaks,'Visible','on');
        selectPeaksButtons{i}=text('Parent',FeatureActionsH,'Units','data','FontSize',10,'VerticalAlignment','Bottom','HorizontalAlignment','Left','String','\bfP','Position',[9.1,nFeatures-i+.15],'BackGroundColor',guiColor,'Tag', int2str(i),'ButtonDownFcn', @selectPeaks,'Visible','on');
        peakH{i}=patch( ...
            'Parent',aH, ...
            'XData',[0,0,0,0], ...
            'YData',[0,0,0,0], ...
            'LineWidth',1, ...
            'FaceColor',peakColor, ...
            'EdgeColor',peakColor, ...
            'FaceAlpha',peakAlphaDefault, ...
            'EdgeAlpha',peakAlphaDefault*2);
        
    end
end

%set(mainFeatureButtons,'SelectedObject',[]);
set(PPSbuttons,'SelectedObject',[]);
%set(FeatureStatsButtons,'SelectedObject',[]);
PPSstate='';
selectSong();
uistack(fig,'top');
%%

%PPS ACTIONS
    function selectPPS(hObject,eventdata)
        PPSstate=get(eventdata.NewValue,'Tag'); % Get Tag of the selected object.
        switch PPSstate
            case 'play'
                % Code for when radiobutton1 is selected.
                playPlayer()
            case 'pause'
                % Code for when radiobutton2 is selected.?
                pausePlayer()
            case 'stop'
                % Code for when togglebutton1 is selected.
                stopPlayer()
                
            otherwise
                % Code for when there is no match.
        end
    end
%ZOOM
    function zoomAxes(hObject,eventdata)
        if not(ishandle(fig))
            return
        end
        
        if strcmp(get(hObject,'Tag'),'in')
            zoomFactor=zoomFactorDefault;
        else
            zoomFactor=1/zoomFactorDefault;
        end
        sliderHLim=get(sliderH,'XData');
        
        sliderWidth=zoomFactor*(sliderHLim(2)-sliderHLim(1));
        if get(followPointerButton,'Value') %is activated -> go to pointer location
            sliderPosition=get(pointerH,'XData');
            sliderPosition=mean(sliderPosition(1:2))/(xlim(2)-xlim(1));
        else
            sliderPosition=mean(sliderHLim);
        end
        
        sliderHLim([1,4])=max(0,min(sliderPosition-sliderWidth/2, 1-sliderWidth));
        sliderHLim([2,3])=min(1,sliderHLim(1)+sliderWidth);
        set(sliderH,'XData',sliderHLim);
        set(aH,'Xlim',xlim(1)+sliderHLim(1:2)*(xlim(2)-xlim(1)));
        
    end

%SLIDE
    function slideAxes(hObject,eventdata)
        %used while mouse button is held down
        CurrentPoint= get(sliderAxes, 'CurrentPoint');
        axesLim=get(aH,'Xlim');
        
        if all(CurrentPoint(1,1:2)<=1) && all(CurrentPoint(1,1:2)>=0) %mouse moved inside sliderAxes
            
            
            sliderWidth=(axesLim(2)-axesLim(1))/(xlim(2)-xlim(1));
            
            sliderLim(1)=max(0, min(CurrentPoint(1,1)-pointOnSlider, 1-sliderWidth));
            sliderLim(2)=min(1,sliderLim(1)+sliderWidth);
            set(sliderH,'XData',[sliderLim(1), sliderLim(2), sliderLim(2), sliderLim(1)]);
            
            set(aH,'Xlim',xlim(1)+sliderLim*(xlim(2)-xlim(1))); % slider at the center of mouse movement
            
        end
    end

    function stopSlide(hObject,eventdata)
        set(fig,'WindowButtonMotionFcn', '');
        
        if strcmp(PPSstate,'play') && get(followPointerButton,'Value')
            return
        end
        
        CurrentPoint= get(sliderAxes, 'CurrentPoint');
        axesLim=get(aH,'Xlim');
        
        if all(CurrentPoint(1,1:2)<=1) && all(CurrentPoint(1,1:2)>=0) %mouse button released inside sliderAxes
            
            
            sliderWidth=(axesLim(2)-axesLim(1))/(xlim(2)-xlim(1));
            
            sliderLim(1)=max(0, min(CurrentPoint(1,1)-pointOnSlider, 1-sliderWidth));
            sliderLim(2)=min(1,sliderLim(1)+sliderWidth);
            set(sliderH,'XData',[sliderLim(1), sliderLim(2), sliderLim(2), sliderLim(1)]);
            
            set(aH,'Xlim',xlim(1)+sliderLim*(xlim(2)-xlim(1))); % slider at the center of mouse movement
            
        end
    end

%PLAY
    function playPlayer(hObject, eventdata)
        
        if not(ishandle(fig))
            return
        end
        if not(isplaying(player))
            
            %start playing from the last frame beginning before CurrentSample (with novelty curves this is the frame center) position
            pointerAlpha=pointerAlphaDefault;
            set(pointerH,'FaceAlpha',pointerAlpha,'EdgeAlpha',pointerAlpha);
            set(smallPointerH,'XData',(CurrentSample-1)/player.TotalSamples*[1,1]);
            play(player,CurrentSample);
        else
            %do nothing
        end
    end

%PAUSE
    function pausePlayer(hObject, eventdata)
        
        if not(ishandle(fig))
            return
        end
        
        %if strcmp(PPSstate,'play')
        pause(player);

        set(pointerH,'XData',framePos([1,2,2,1],CurrentFrame));
        set(smallPointerH,'XData',(CurrentSample-1)/player.TotalSamples*[1,1]);
        drawnow
    end

%STOP
    function stopPlayer(hObject, eventdata)
        
        if not(ishandle(fig))
            return
        end
        stop(player);
        CurrentSample=get(player,'CurrentSample');
        CurrentFrame=1;
        
        set(pointerH,'XData',xlim([1,1,1,1]));
        set(smallPointerH,'XData',[0,0]);
        drawnow
        
    end

    function startDragFcn(varargin)
        
        if not(ishandle(fig))
            return
        end
        
        CurrentPointAxes=get(aH, 'CurrentPoint');
        CurrentPointSlider=get(sliderAxes,'CurrentPoint');
        
        if CurrentPointAxes(1,1)>=xlim(1) && CurrentPointAxes(1,1)<=xlim(2) && CurrentPointAxes(1,2)>=ylim(1) && CurrentPointAxes(1,2)<=ylim(2) %mouse clicked inside aH
            if strcmp(PPSstate,'play')
                pause(player);
            end
            
            
            
            
        [tmp, CurrentFrame]=min(framePos(1,:)<(xlim(1)+CurrentPointAxes(1,1)));
            
            %CurrentFrame=length(find(framePos(1,:)<=CurrentPointAxes(1,1)));
            %pointerX=[framePos(1,CurrentFrame),framePos(2,CurrentFrame),framePos(2,CurrentFrame),framePos(1,CurrentFrame)];
            pointerAlpha=pointerAlphaDefault;
            set(pointerH,'XData',framePos([1,2,2,1],CurrentFrame),'FaceAlpha',pointerAlpha,'EdgeAlpha',pointerAlpha);
            set(fig, 'WindowButtonMotionFcn', @draggingFcn)
            set(fig, 'WindowButtonUpFcn', @stopDragFcn)
        elseif all(CurrentPointSlider(1,1:2)<=1) && all(CurrentPointSlider(1,1:2)>=0) %mouse clicked inside sliderAxes
            %disable sliding during playback if followPointerButton
            %selected
            if strcmp(PPSstate,'play') && get(followPointerButton,'Value')
                return
            else
                axesLim=get(aH,'Xlim');
                sliderWidth=(axesLim(2)-axesLim(1))/(xlim(2)-xlim(1));
                
                sliderLim=get(sliderH,'XData');
                sliderLim=sliderLim(1:2); %get slider position
                
                if sum(sign(CurrentPointSlider(1,1)-sliderLim))~=0  %if the current point is not inside the slider -> move the slider to the desider position
                    
                    sliderLim(1)=max(0, min(CurrentPointSlider(1,1)-sliderWidth/2, 1-sliderWidth));
                    sliderLim(2)=min(1,sliderLim(1)+sliderWidth);
                    set(sliderH,'XData',[sliderLim(1), sliderLim(2), sliderLim(2), sliderLim(1)]);
                    
                    set(aH,'Xlim',xlim(1)+sliderLim*(xlim(2)-xlim(1))); % slider at the center of mouse click
                    set(smallPointerH,'XData',(player.CurrentSample-1)/player.TotalSamples*[1,1]);
                end
                pointOnSlider=CurrentPointSlider(1,1)-sliderLim(1);
                set(fig, 'WindowButtonMotionFcn', @slideAxes);
                set(fig, 'WindowButtonUpFcn', @stopSlide)
            end
            
        else
            return
        end
        
        
        
    end

    function draggingFcn(varargin)
        %used while mouse button is held down
        CurrentPoint= get(aH, 'CurrentPoint');
        
        if CurrentPoint(1,1)<xlim(1) || CurrentPoint(1,1)>xlim(2) || CurrentPoint(1,2)<ylim(1) || CurrentPoint(1,2)>ylim(2)
            return
        end
        
        
        [tmp, CurrentFrame]=min(framePos(1,:)<(xlim(1)+CurrentPoint(1,1)));
        %pointerX =framePos([1,2,2,1],CurrentFrame);% [framePos(1,CurrentFrame), framePos(2,CurrentFrame), framePos(2,CurrentFrame), framePos(1,CurrentFrame)];
        
        %CurrentFrame=length(find(framePos(1,:)<=CurrentPoint(1,1)));
        
        %pointerX=[framePos(1,CurrentFrame),framePos(2,CurrentFrame),framePos(2,CurrentFrame),framePos(1,CurrentFrame)];
        set(pointerH,'XData',framePos([1,2,2,1],CurrentFrame));
        set(smallPointerH,'XData',CurrentPoint(:,1));
    end

    function stopDragFcn(varargin)
        set(fig,'WindowButtonMotionFcn', '');
        
        %pointerX=[framePos(1,CurrentFrame),framePos(2,CurrentFrame),framePos(2,CurrentFrame),framePos(1,CurrentFrame)];
        set(pointerH,'XData',framePos([1,2,2,1],CurrentFrame));
        set(smallPointerH,'XData',CurrentFrame/length(framePos)*[1,1]);
        
        %get current sample from current frame
        CurrentSample=max(1,floor(Fs*(framePos(1,CurrentFrame)-xlim(1))));
        %%if play button is activated
        if strcmp(PPSstate,'play')
            play(player,CurrentSample);
        end
        
    end

    function selectSong(varargin)
        % select a song from miraudio struct and plot it in aH. Update also
        % possible curve and peak data.
        
        if not(ishandle(fig))
            return
        end
        
        if isequal(class(player),'audioplayer')
            stopPlayer();  % No selection
            player=[];
            java.lang.Runtime.getRuntime.gc; %Java garbage collection
            
        else
            CurrentSample=1;
        end
        
        if isequal(get(colorPatch,'Visible'),'on')
            openedFeature=str2num(get(colorPatch,'Tag'));
            set(colorPatch,'Visible','off');
            set(featureText{openedFeature},'Visible','on');
            set(selectFeatureButtons{openedFeature},'XData',[.1,.9,.1],'YData',nFeatures-openedFeature+[.15,.5,.85],'ButtonDownFcn',@selectFeature);
            set(selectPeaksButtons{openedFeature},'ButtonDownFcn',@selectPeaks);
        end
        
        songInd=get(audioPopupmenuH, 'Value');
        song = miraudio(features.songNames{songInd});    
        
        %start and end in seconds
        xlim=get(song,'Pos');
        xlim = xlim{1}{1}([1 end]);
        ylim=[0,1];
        Fs=get(song,'Sampling');
        Fs = Fs{1};
        
        try
            player = audioplayer(mirgetdata(song), Fs);
        catch exception
            fixException(exception)
        end
        
        %player.StartFcn = 'startPlaying';
        player.TimerPeriod = maxFrameUpdateFrequency;
        player.TimerFcn = 'whilePlaying';
        
        
        set(aH,'Xlim',xlim);
        set(sliderAxes,'Xlim',[0,1]);
        
        songImage=mirgetdata(miraudio(song,'Sampling',downSampleRate));
        %songImage=.5+.5*(songImage-mean(songImage))/max(abs(songImage));
        
        songImage=(songImage-min(songImage))/(max(songImage)-min(songImage))-.5; %normalize to [-.5,.5], mean=0
        
        songImagePos=(0:length(songImage)-1)./downSampleRate;
        
        %could handle stereo wave -> take mean across channels?
        
        set(songThumbnailH,'XData',(0:length(songImage)-1)/(length(songImage)-1),'YData', songImage+.5, 'Color', songColor);
        set(songH,'XData',songImagePos,'YData',.75*(ylim(2)-ylim(1))*songImage+mean(ylim),'Color',songColor);
        selectedFeatureInds=find(selectedFeatures~=2);  %2 means hidden feature
        
        if isempty(selectedFeatureInds)
            CurrentSample=1;
            setAudioBackground();
            
        else
            for sFI=1:length(selectedFeatureInds)
                selectFeature([],[],selectedFeatureInds(sFI),selectedFeatures(selectedFeatureInds(sFI)));
            end
            
        end
        
        
        
    end


    function selectFeature(hObject, eventdata, featureInd, featureState)
        songInd=get(audioPopupmenuH, 'Value');
        
        if isequal(get(colorPatch,'Visible'),'on')
            openedFeature=str2num(get(colorPatch,'Tag'));
            set(colorPatch,'Visible','off');
            set(featureText{openedFeature},'Visible','on');
            set(selectFeatureButtons{openedFeature},'XData',[.1,.9,.1],'YData',nFeatures-openedFeature+[.15,.5,.85],'ButtonDownFcn',@selectFeature);
            set(selectPeaksButtons{openedFeature},'ButtonDownFcn',@selectPeaks);
        end
        
        set(noDataText,'Visible', 'off');
        
        if isempty(hObject)
            if isempty(features.data{featureInd}{songInd}{1})
                set(noDataText,'Visible', 'on');
                set(featureH{featureInd},'XData',[0,0],'YData',[0,0],'Visible','off');
                if ~mainFeature
                    setAudioBackground();
                end
                return
            else
                setFeature([],[],featureInd,featureState);
            end
            
        else
            
            featureInd=str2num(get(hObject,'Tag'));
            
            if features.emptysong(songInd)
                showFeatureStats([],[],featureInd);
                if ~mainFeature
                    setAudioBackground();
                end
                return
                
            end
            
            
            set(colorPatch,'Visible','on','Vertices',[colorVertices(:,1),colorVertices(:,2)+nFeatures-featureInd],'Faces',colorFaces,'FaceVertexCdata',colorCdata,'Tag',get(hObject,'Tag'),'ButtonDownFcn',@setFeature);
            set(selectFeatureButtons{featureInd},'XData',[.1,.9,.9],'YData',nFeatures-featureInd+[.5,.15,.85],'ButtonDownFcn',@setFeature);
            set(featureText{featureInd},'Visible','off');
            
            showFeatureStats([],[],featureInd);
            
            set(fig,'WindowButtonDownFcn','','WindowButtonUpFcn', '');
        end
    end



    function setFeature(hObject, eventdata, featureInd,featureState)
        %if isempty(hObject)
        %    featureColor=featureColors(featureState,:);
        if not(isempty(hObject))
            clickedPoint=get(FeatureActionsH,'CurrentPoint');
            featureInd=str2num(get(hObject,'Tag'));
            
            clickedPoint=ceil(clickedPoint(1,1));
            
            featureStates=[0,1,1,2,2,3,4,5,6];
            featureState=featureStates(clickedPoint);
            
            set(colorPatch,'Visible','off');
            set(featureText{featureInd},'Visible','on');
            set(selectFeatureButtons{featureInd},'XData',[.1,.9,.1],'YData',nFeatures-featureInd+[.15, .5,.85],'ButtonDownFcn',@selectFeature);
        end
        
        if featureState>0
            set(selectFeatureButtons{featureInd},'FaceColor',featureColors(featureState,:));
        end
        
        if featureState==0
            set(selectFeatureButtons{featureInd},'ButtonDownFcn',@selectFeature);
            %only hide the colorpatch etc.
        elseif featureState==1
            selectedFeatures(featureInd)=1;
            if mainFeature>0 && mainFeature~=featureInd  %Hide the previous main feature. There can be only one at a time.
                set(selectFeatureButtons{mainFeature},'FaceColor',featureColors(2,:));
                selectedFeatures(mainFeature)=1;
                %feature{mainFeature}=[];
                %featureCurvePos{mainFeature}=[];
                set(featureH{mainFeature},'XData',[0,0],'YData',[0,0],'Visible','off');
            end
            mainFeature=featureInd;
            setMainFeature(featureInd);
        elseif featureState==2
            selectedFeatures(featureInd)=2;
            %feature{featureInd}=[];
            %featureCurvePos{featureInd}=[];
            set(featureH{featureInd},'XData',[0,0],'YData',[0,0],'Visible','off');
                        if featureInd==mainFeature
                mainFeature=0;
            end
        elseif featureState>=3
            selectedFeatures(featureInd)=featureState;
            setSecondaryFeature(featureInd,featureState);
            if featureInd==mainFeature
                mainFeature=0;
            end
            
        end
        %set(FeatureActionsH,'ButtonDownFcn','');
        set(fig,'WindowButtonDownFcn', @startDragFcn);
    end


    function setSecondaryFeature(featureInd,featureState)
        if not(ishandle(fig))
            return
        end
        
        songInd=get(audioPopupmenuH, 'Value');
        %[feature{featureInd} featureCurvePos{featureInd}, featureSongDistribution]= ...
        %    getSongFeatureData(f,featureFields{featureInd},songInd,featureValueRange(featureInd,:),nBins,smoothingFactor);
        %showFeatureStats([],[],featureInd, featureSongDistribution);
        
        %center feature to [-.5,.5]
        if ~features.isSongLevel(featureInd) && features.minsong(featureInd,songInd) ~= features.maxsong(featureInd,songInd)
            %minValue=features.minsong(featureInd,songInd);
            %maxValue=features.maxsong(featureInd,songInd);
            
            set(featureH{featureInd},...'XData',mean(features.framePos{featureInd}{songInd}{1}), ...
                ...'YData',((features.data{featureInd}{songInd}{1}-minValue)/(maxValue-minValue)-.5)*.9*(ylim(2)-ylim(1))+mean(ylim), ...
                'Color',featureColors(featureState,:),'LineWidth',1,'Visible','on');            
        else
            

            set(featureH{featureInd},'XData',xlim,'YData',mean(ylim)*ones(1,2),'Color',featureColors(featureState,:),'LineWidth',1,'Visible','on');
        end
        
    end


    function setMainFeature(featureInd)
        %Select a feature in popupmenu. Plot the feature curve. At this
        %point only mirscalars are allowed.
        if not(ishandle(fig))
            return
        end
        
        
        %%if play button is activated
        if strcmp(PPSstate,'play')
            pausePlayer();
        end
        
        songInd=get(audioPopupmenuH, 'Value');
        showFeatureStats([],[],featureInd);
        
        
        %compute feature information for visualization
        
        
        if ~features.isSongLevel(featureInd) && features.minsong(featureInd,songInd) ~= features.maxsong(featureInd,songInd)
            framePos=features.framePos{featureInd}{songInd}{1};
            %add 5% overhead in aH
            ylim=[features.minsong(featureInd,songInd),features.maxsong(featureInd,songInd)];
            ylim=ylim+.05*[ylim(1)-ylim(2),ylim(2)-ylim(1)];
        else
            framediff=maxFrameUpdateFrequency; %play audio normally
            framePos=[xlim(1):framediff:xlim(2);(xlim(1)+framediff):framediff:(xlim(2)+framediff)];
            %feature will be in the center of aH (in y-direction)
            ylim=[min([0,features.minsong(featureInd,songInd)]),2*max([0,features.maxsong(featureInd,songInd)])];
            if all(ylim)==0
                ylim=ylim+[-.05,.05];
            end
        end
        
        CurrentSample=get(player,'CurrentSample');
        [tmp, CurrentFrame]=min(framePos(1,:)<(xlim(1)+CurrentSample/Fs));
        
        set(aH,'Ylim',ylim);
        
        if ~features.isSongLevel(featureInd) && features.minsong(featureInd,songInd) ~= features.maxsong(featureInd,songInd)
            if iscell(features.data{featureInd}{songInd}{1})
                xdata = [];
                ydata = [];
                for k = 1:length(features.data{featureInd}{songInd}{1})
                    xdata = [xdata repmat(mean(framePos(:,k)),[1 length(features.data{featureInd}{songInd}{1}{k})])];
                    ydata = [ydata features.data{featureInd}{songInd}{1}{k}(:)'];
                end
                linestyle = 'none';
                marker = '+';
            else
                xdata = mean(framePos);
                ydata = features.data{featureInd}{songInd}{1};
                linestyle = '-';
                marker = 'none';
            end
            set(featureH{featureInd},'XData', xdata, ...
                'YData',ydata,'Color','k','LineWidth',1,'Visible','on','LineStyle',linestyle,'Marker',marker);
        else

            set(featureH{featureInd},'XData', xlim, ...
                'YData',mean(ylim) * ones(1,2),'Color','k','LineWidth',1,'Visible','on');            
        end
        

        
        %scale secondary features and peaks to ylim
        for j=1:length(selectedFeatures)
            if selectedFeatures(j)>2
                
                minValue=min(features.data{j}{songInd}{1});
                maxValue=max(features.data{j}{songInd}{1});
                
                set(featureH{j},'YData',((features.data{j}{songInd}{1}-minValue)/(maxValue-minValue)-.5)*.9*(ylim(2)-ylim(1))+mean(ylim));
            end
            if selectedPeaks(j)>1
                set(peakH{j},'YData',[zeros(2,length(features.peaks{j}{songInd}));0.945*(ylim(2)-ylim(1))*repmat(features.peaks{j}{songInd},2,1)] + ylim(1));
            end
        end
        
        uistack(featureH{featureInd},'top');
        
        set(songH,'YData',.75*(ylim(2)-ylim(1))*songImage+mean(ylim));
        
        set(pointerH,'XData',framePos([1,2,2,1],CurrentFrame),'YData',ylim([1,1,2,2])+[-.1,-.1,.1,.1]);
        set(smallPointerH,'XData',(CurrentSample-1)/player.TotalSamples*[1,1]);
        
        %the pointer update freq will be an integer [1,2,3,...]
        %multiplication of the feature framediff, at least
        %maxFrameUpdateFrequency (to reduce the computational burden with short frame lengths)
        %nFramesWhenUpdate=ceil(maxFrameUpdateFrequency/framediff);  
        %updateFreq=framediff*nFramesWhenUpdate;
        %player.TimerPeriod = updateFreq;
        
        %%if play button is activated
        if strcmp(PPSstate,'play')
            playPlayer();
        end
        
    end

    function setAudioBackground()
        framediff=maxFrameUpdateFrequency; %play audio normally
        framePos=[xlim(1):framediff:xlim(2);(xlim(1)+framediff):framediff:(xlim(2)+framediff)];
        
        %pointerX=[framePos(1,CurrentFrame),framePos(2,CurrentFrame),frameP
        %os(2,CurrentFrame),framePos(1,CurrentFrame)];
        set(aH,'Ylim',ylim);
        set(songH,'YData',.75*(ylim(2)-ylim(1))*songImage+mean(ylim));
        
        pointerY=[ylim(1), ylim(1), ylim(2), ylim(2)]+[-.1,-.1,.1,.1];
        set(pointerH,'XData',framePos([1,2,2,1],CurrentFrame),'YData',pointerY);
        set(smallPointerH,'XData',(CurrentSample-1)/player.TotalSamples*[1,1]);
    end


    function selectPeaks(hObject, eventdata, featureInd,peakState)
        %select peaks in popupmenu. Called usually each time when peak popupmenu
        %called or song selected (if there are peaks at all in the feature set)
        
        if not(ishandle(fig))
            return
        end
        
        
        if isempty(hObject)
            setPeaks([],[],featureInd,peakState);
        else
            if isequal(get(colorPatch,'Visible'),'on')
                openedFeature=str2num(get(colorPatch,'Tag'));
                set(colorPatch,'Visible','off');
                set(featureText{openedFeature},'Visible','on');
                set(selectFeatureButtons{openedFeature},'XData',[.1,.9,.1],'YData',nFeatures-openedFeature+[.15,.5,.85],'ButtonDownFcn',@selectFeature);
                set(selectPeaksButtons{openedFeature},'ButtonDownFcn',@selectPeaks);
            end
            peakInd=str2num(get(hObject,'Tag'));
            set(colorPatch,'Visible','on','Vertices',[peakColorVertices(:,1),peakColorVertices(:,2)+nFeatures-peakInd], ...
                'Faces',peakColorFaces,'FaceVertexCdata',peakColorCdata,'Tag',get(hObject,'Tag'),'ButtonDownFcn',@setPeaks);
            %set(FeatureActionsH,'ButtonDownFcn',@setFeature);
            set(selectPeaksButtons{peakInd},'ButtonDownFcn',@setPeaks);
            set(featureText{peakInd},'Visible','off');
            songInd=get(audioPopupmenuH, 'Value');
            %[tmp,tmp, featureSongDistribution]=getSongFeatureData(f,featureFields{featureInd},songInd,featureValueRange(featureInd,:),nBins);
            %showFeatureStats([],[],featureInd, featureSongDistribution);
            
            set(fig,'WindowButtonDownFcn','','WindowButtonUpFcn', '');
        end
        
    end

    function setPeaks(hObject, eventdata, featureInd,peakState)
        %if isempty(hObject)
        %    peakColor=featureColors(peakState,:);
        %else
        if not(isempty(hObject))
            clickedPoint=get(FeatureActionsH,'CurrentPoint');
            featureInd=str2num(get(hObject,'Tag'));
            
            clickedPoint=ceil(clickedPoint(1,1));
            
            peakStates=[0,0,0,2,3,4,5,1,1,0];
            peakState=peakStates(clickedPoint);
            
            set(colorPatch,'Visible','off');
            set(featureText{featureInd},'Visible','on');
            set(selectPeaksButtons{featureInd},'ButtonDownFcn',@selectPeaks);
            
            if peakState>0
                set(selectPeaksButtons{featureInd},'BackGroundColor',peakColors(peakState,:));
            end
            
        end
        
        if peakState==0
            set(selectFeatureButtons{featureInd},'ButtonDownFcn',@selectPeaks);
            %only hide the colorpatch etc.
        elseif peakState==1
            selectedPeaks(featureInd)=1;
            set(peakH{featureInd},'XData',[0,0],'YData',[0,0],'Visible','off');
        elseif peakState>=2
            songInd=get(audioPopupmenuH, 'Value');
            selectedPeaks(featureInd)=peakState;
            %update patch data: each column in matrices are one patch
            peakPos=features.framePos{featureInd}{songInd}{1}(:,features.peakPos{featureInd}{songInd});

            if strcmp(PPSstate,'play')
                pausePlayer();
            end
            set(peakH{featureInd},'XData',peakPos([1,2,2,1],:),'YData',[zeros(2,size(peakPos,2))+ylim(1);.9*(ylim(2)-ylim(1))*repmat(features.peaks{featureInd}{songInd},2,1) + abs(min(features.data{featureInd}{songInd}{1})-ylim(1))],...%min(features.data{featureInd}{songInd}{1})], ...
                'FaceColor',peakColors(peakState,:),'EdgeColor',peakColors(peakState,:),'Visible','on');
            
            %%if play button is activated
            if strcmp(PPSstate,'play')
                playPlayer();
            end
            
            
        end
        %set(FeatureActionsH,'ButtonDownFcn','');
        set(fig,'WindowButtonDownFcn', @startDragFcn, 'WindowButtonUpFcn', @stopDragFcn);
        
        
    end



    function showFeatureStats(hObject,eventdata,featureInd)
        songInd=get(audioPopupmenuH, 'Value');
        
        if not(isempty(hObject))
            featureInd=str2num(get(eventdata.NewValue, 'Tag'));
            %[~,~,featureSongDistribution]=getSongFeatureData(f,featureFields{featureInd},songInd,featureValueRange(featureInd,:),nBins,smoothingFactor);
        end
        
        ticklabels{1}=num2str(features.valueRange(featureInd,1));
        ticklabels{2}=num2str(mean(features.valueRange(featureInd,:)));
        ticklabels{3}=num2str(features.valueRange(featureInd,2));
        
        if length(ticklabels{1})>4
            ticklabels{1}=num2str(features.valueRange(featureInd,1),'%1.2e');
        end
        if length(ticklabels{2})>4
            ticklabels{2}=num2str(mean(features.valueRange(featureInd,:)),'%1.2e');
        end
        if length(ticklabels{3})>4
            ticklabels{3}=num2str(features.valueRange(featureInd,2),'%1.2e');
        end
        
        if features.emptysong(songInd)
            set(noDataText,'Visible','on')
        else
            set(noDataText,'Visible','off')
        end
        set(DistPanel,'Title',upper(features.names{featureInd}));
        set(distAxes,'Xtick',[0,.5,1],'XTickLabel',ticklabels);
        set(featureDistPatch,'YData',[features.distribution(featureInd,:),0,0]);
        if ~isempty(features.songDistributions{featureInd})
            set(songDistPatch,'YData',[features.songDistributions{featureInd}(songInd,:),0,0]/2);
        end
        
    end

    function fixException(exception)
        message=exception.getReport('basic');
        expectedMsg='java.lang.OutOfMemoryError: Java heap space';
        if strcmp(message(end+1-length(expectedMsg):end),expectedMsg)
            %Increase heap size. See tutorial at
            %http://www.mathworks.com/support/solutions/en/data/1-18I2C/
            heaps=2.^(1:20);
            currentHeap=java.lang.Runtime.getRuntime.maxMemory/1000000;
            [tmp,heapInd]=min(dist(currentHeap,heaps));
            
            if heapInd==size(heaps,2)
                disp('Too much heap space already. Aborting...');
                return
            end
            
            currentHeap=heaps(heapInd);
            
            %These will be used below
            totalMemory=heaps(6);
            maxMemory=max(heaps(heapInd+1),heaps(9));
            
            question=sprintf('Maximum Java heap space in MATLAB is currently too little for playback (%dMB). Increase the maximum Java heap size to %dMB (requires restarting MATLAB)?', ...
                currentHeap, heaps(9));
            
            heapButton = questdlg(question, ...
                'Java Heap Space','Yes','No','No');
            switch heapButton
                case 'Yes',
                    up=userpath;
                    if up(end)==':'
                        up(end)=[];
                    end
                    
                    try
                        javaver=version('-java');
                        javaver=strread(javaver,'%s','delimiter','_');
                        javaver=strread(javaver{1},'%s','delimiter',' ');
                        javaver=strread(javaver{2},'%d','delimiter','.');
                    catch javaVerException
                        javaver=strread(strtrim(input('Could not parse Java version. Give the correct java version number in the form similar to 1.6.0 ( type version(''-java'') to see your version).'),'s'),'%d','delimiter','.');
                    end
                    
                    if javaver(1)<1 || javaver(1)==1 && javaver(2)<2 || javaver(1)==1 && javaver(2)==2 && javaver(3)<2 % java version must be 1.2.2 or later
                        display('Java version must be 1.2.2 or later. Aborting...');
                        return
                    end
                    
                    
                    if all(javaver==[1;1;8]) %different strings
                        heapbeginning='-';
                    else
                        heapbeginning='-X';
                    end
                    
                    javaopts=fopen(strcat(up,'/java.opts'),'w');
                    
                    
                    disp('Writing java.opts file to be loaded in MATLAB startup (see the userpath).');
                    %Quartz is faster than the rendering pipeline is provided by
                    %Sun
                    fprintf(javaopts,'-Dapple.awt.graphics.UseQuartz=true\n%sms%dm\n%smx%dm',heapbeginning,totalMemory,heapbeginning,maxMemory);
                    fclose(javaopts);
                    quit
                    
                case 'No',
                    disp('Aborting...');
                    return
            end
            
            
        else
            throw(exception);
        end
    end

end
