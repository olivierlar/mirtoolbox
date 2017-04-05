function mirplayer(FEATURES,SELECT)
%MIRPLAYER is a graphical player interface MIRtoolbox objects
%   MIRPLAYER(FEATURES) displays the interface. FEATURES is a mirtoolbox object
%   such as mirstruct or mirscalar, or a nested MATLAB structure of mirtoolbox
%   objects.
%
%   MIRPLAYER(FEATURES,SELECT) displays only songs specified by an array 
%   of song indices SELECT.
%
%   Examples:
%
%   The examples assume that the folder MIRtoolboxX.X.X/MIRToolboxDemos
%   is in the MATLAB path.
%   
%   Display a feature
%   FEATURES = mirbrightness('ragtime','Frame');
%   MIRPLAYER(FEATURES)
%
%   Display a multiple features
%   clear FEATURES;
%   FEATURES.brightness = mirbrightness('ragtime','Frame');
%   FEATURES.filterbank = mirfilterbank('ragtime');
%   MIRPLAYER(FEATURES)
%
%   or
%   FEATURES = mirfeatures('ragtime');
%   MIRPLAYER(f)
%
%   Display audio file
%   a = miraudio('ragtime');
%   MIRPLAYER(a)
%
%   Display feature(s) computed for multiple audio files
%   cd MIRToolboxDemos/; % specify here any folder that contains audio files
%   FEATURES = mirbrightness('Folder','Frame');
%   MIRPLAYER(FEATURES)
%
%   or display only the first and the fifth audio file
%   SELECT = [1,5];
%   MIRPLAYER(FEATURES,SELECT)
%
%Certain limitations in the functionality:
%- Not optimized for displaying song level features (such as those 
%   extracted without the 'Frame' parameter).
%
%   See also MIRPLAY.


if ischar(FEATURES) %session data given in a file
    if (exist(strcat(FEATURES,'/featureInfo.mat'),'file') ~=2) || ~exist(strcat(FEATURES,'/1.wav'),'file')
        error('Session data file not found in a MATLAB file.');
    end
    disp('Loading the session. Please wait...');
    load(strcat(FEATURES,'/featureInfo.mat'));
    if ~exist('feat_info','var') || ~exist('songNames','var') || ~exist('xlims','var') || ~exist('songSampling','var')
        error('Session feature data file is not compatible with MIRplayer');
    end
    nBins=size(feat_info.distribution,2);
    songs=1:length(songNames);
else
    %the first argument should be the feature set
    if ~isstruct(FEATURES) && ~iscell(FEATURES) && ~isa(FEATURES,'mirdata')
        error('The first input argument should be struct or mirdata variable.');
    end
    
    if nargin==2
        songs=int8(unique(SELECT));
        
        %if any(songs)<1 || any(songs)>length(songNames)
        %    error('The second input argument should be an integer array of songs to be included in the analysis.');
        %end
    else
        songs=[];
    end
    
    nBins=100; %resolution of the feature distributions (for visualization)
    smoothingFactor=2;%round(max(3,nBins/10)); %for median filtering the distribution
    feat_info=getFeatureInfo(FEATURES, nBins,smoothingFactor,songs);
    
    %TODO: Audio could also be shown without features. Just to be able to check
    %out the songs before extracting features.
    if isempty(feat_info)
        error('Please provide a feature set with at least one mirscalar type of feature.');
    end
    
end


global fig
global Fs
global pointerH
global playheads
global player
global CurrentSample
global CurrentFrame
global framePos
global xlim
global PlayPauseButton
global playIcon
global PPSstate
global smallPointerH
global sliderH
global aH
global followPointerButton
global featureAxes
global frameSummary
global loopingButton

%framediff=0;
song=0;
songImage=0;
playheads=[];
featureAxes={};
selectedFeatures=[];
%pointerY=[];
ylim=[0,1];
%CurrentFrameStartPosition=0;

pointOnSlider=[];

%mirchunklim=5000000;

nFeatures=length(feat_info.names);

guiColor=[.93,.93,.93];
pointerColor1='k'; %Black
%peakColor='r'; %Red
%featureColors=[0,0,0;guiColor;0,.6,0;.9,.7,0;0,0,.9;.8,0,0];%0,.5,0;.9,.6,0;0,0,.8;.6,0,0]; %available rgb values for plotting the features
%peakColors=[guiColor;0,.5,0;.9,.6,0;0,0,.8;.6,0,0]; %available rgb values for plotting the features
pointerColor=pointerColor1;
playheadAlpha=.3;
%peakAlphaDefault=.1;
%downsample a to save memory when plotting audio
%downSampleRate=1000; %Could be nice resolution

songColor=[.85,.85,.85];
maxFrameUpdateFrequency=.05; %seconds
zoomFactorDefault=2/3;
%feature=cell(nFeatures,1);
%featureCurvePos=cell(nFeatures,1);
%selectFeatureButtons=cell(nFeatures,1);
%selectPeaksButtons=cell(nFeatures,1);
%mainFeature=0;
%selectedFeatures=2*ones(nFeatures,1);
%selectedFeatureInd=0;

%selectedPeaks=ones(nFeatures,1);

frameSummary=[];

%create play, pause and stop icons for buttons
%s=30;
mp = get(0, 'MonitorPositions');
s=mp(1,:);%get(0,'ScreenSize');
s=round(s(3)/150)*2; %size of the rectangular area
stopIcon=zeros(s,s,3);
pauseIcon=zeros(s,s,3);
pauseIcon(:,s/2-1:s/2+1,:)=NaN;
tri=triu(ones(s/2,s/2),1);
playIcon=zeros(s,s); playIcon(1:s/2,1:2:s)=tri;
playIcon(1:s/2,2:2:s)=tri;
playIcon=repmat(playIcon+playIcon(s:-1:1,:),[1,1,3]);
playIcon(playIcon==1)=NaN;

%read logo
try
    S=dbstack('-completenames');
    logo=double(imread(regexprep(S(1).file,[S(1).name,'.m'],'../UsersManual/logo.png'), 'PNG'));
    %logo=squeeze(logo(:,:,1));
    %logo = imresize(logo, .5,'bilinear');
    logo(logo~=255)=.87;%guiColor/2;
    logo(logo==255)=guiColor(1);
    
    
    
catch
    %didn't find the logo file
end


%% CREATE GUI
%'Units', 'normalized', ...
fig    =   figure(...       % the main GUI figure
    'MenuBar','fig', ...
    'Toolbar','none', ...
    'HandleVisibility','callback', ...    
    'OuterPosition', mp(1,:), ...
    'PaperUnits', 'inches', ...
    'PaperOrientation', 'landscape', ...
    'PaperSize', [16,9], ...
    'Name', mfilename, ...
    'NumberTitle','off', ...
    'WindowButtonDownFcn', @startDragFcn, ...
    'Color', guiColor);
%menuItems = uimenu('Parent',fig,'Label','File');
%uimenu(menuItems,'Label','Print','Callback','printpreview');
%uimenu(menuItems,'Label','Save session','Callback',@saveSession);
%uimenu(menuItems,'Label','Quit session','Separator','on','Accelerator','Q','Callback',@quitSession);

%     function saveSession(hObject,eventdata)
%         folderName=inputdlg('Folder name','Save session');
%         if isempty(folderName{1})
%             return
%         else
%             mkdir(folderName{1});
%             save(strcat(folderName{1},'/featureInfo'), 'feat_info','songNames','xlims','songSampling');
%         end
%         
%         for i=1:length(songs)
%             ind=songs(i);
%         end
%         
%         
%     end

%     function quitSession(hObject,eventdata)
%         quitting = questdlg('Save session before quitting?','Quit','Don''t save','Cancel','Save', 'Save');
%         if isempty(quitting) || isequal(quitting,'Cancel')
%             return;
%         elseif isequal(quitting,'Save')
%             java.lang.Runtime.getRuntime.gc; %Java garbage collection
%             saveSession;
%             if ishandle(fig)
%                 close(fig);
%             end
%         else
%             if ishandle(fig)
%                 close(fig);
%             end
%         end
%     end

mainPanelPos=[.16,.0,.8,.8];
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
    'Clipping', 'off', ...
    'HandleVisibility', 'callback', ...
    'Position',[.16,.82,.8,.16], ...
    'BorderType','none', ...
    'BackGroundColor', guiColor, ...
    'Visible','on');
FeaturePanel = uipanel(...
    'Parent', fig, ...
    'Title', 'SELECT FEATURES', ...
    'Units', 'normalized', ...
    'Clipping', 'on', ...
    'HandleVisibility', 'callback', ...
    'Position',[.01, 0.1, mainPanelPos(1)-.05, 0.9], ...
    'BorderType','none', ...
    'BackGroundColor', guiColor, ...
    'Visible','on');
% DistPanel = uipanel(...
%     'Parent', fig, ...
%     'Title', 'FEATURE DISTRIBUTION', ...
%     'FontSize', 8, ...
%     'Units', 'normalized', ...
%     'Clipping', 'off', ...
%     'HandleVisibility', 'callback', ...
%     'Position',[.01,.08,mainPanelPos(1)-.05,.14], ...
%     'BackGroundColor', guiColor, ...
%     'BorderType','none', ...
%     'Visible','on');
outerPos=get(fig,'OuterPosition');
ratioPos=outerPos(3)*mainPanelPos(3)/(outerPos(4)*mainPanelPos(4));
aH      =   axes(...         % the axes for plotting
    'Parent', MainPanel, ...
    'Units', 'normalized', ...
    'HandleVisibility','callback', ...
    'Visible','off', ...
    'Clipping', 'on', ...    
    'Position',[0.05 0.1 0.95 0.9]);
xlabel(aH,'Time (s)');
ylabel(aH,'Feature value');
PPSsize=[.15,.3];

PlayPauseButton  =   uicontrol(...
    'Parent', ControlPanel, ...
    'Style','PushButton', ...
    'CData',playIcon, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0,0.6,0.04,0.4], ...%[0,0.4,.05 .5], ...%[.1,0.85,buttonSize],...
    'Tag','play', ...
    'CallBack',@playPausePlayer);
StopButton  =   uicontrol(...
    'Parent', ControlPanel, ...
    'Style','PushButton', ...
    'CData',stopIcon, ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[0.04,0.6,0.04,0.4], ...%[.06,.4,.05 .5],...
    'Tag','stop', ...
    'CallBack',@stopPlayer);
audioPopupmenuH=   uicontrol(...    % list of available audio
    'Parent', ControlPanel, ...
    'Units','normalized',...
    'Position',[.25 .7 .5 .1],...
    'Callback', @selectSong, ...
    'HandleVisibility','callback', ...
    'String',feat_info.songNames,...
    'TooltipString','Available audio files', ...
    'Style','popupmenu');



pointerH=patch( ...
    'Parent',aH, ...
    'XData',[0,0,0,0], ...
    'YData',[0,0,0,0], ...
    'LineWidth',1, ...
    'FaceColor',pointerColor, ...
    'EdgeColor',pointerColor, ...
    'FaceAlpha',playheadAlpha, ...
    'EdgeAlpha',playheadAlpha);

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
    'Position',[.944 .55 .027 .45], ...%[.1,0.85,buttonSize],...
    'CData',zoomInIcon, ...
    'Tag', 'in', ...
    'TooltipString','Horizontally zoom in the feature axes', ...
    'CallBack',@zoomAxes);
zoomOutButton  =   uicontrol(...
    'Parent', zoomButtons, ...
    'Style','PushButton', ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'Position',[.974 .55 .027 .45], ...%[.1,0.85,buttonSize],...
    'CData',zoomOutIcon, ...
    'Tag', 'out', ...
    'TooltipString','Horizontally zoom out the feature axes', ...
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
    'Position',[0.05,0,0.95,.5]);%, ...
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
    'String', 'Follow playhead', ...
    'TooltipString','Follow playhead position when playing', ...
    'Position',[.82 .55 .12 .5], ...
    'BackGroundColor',guiColor, ...
    'CallBack', @followPointer);%, ...

loopingButton  =   uicontrol(...
    'Parent', zoomButtons, ...
    'Style','CheckBox', ...
    'Units','normalized',...
    'HandleVisibility','callback', ...
    'String', 'Looping', ...
    'TooltipString','Loop play in the selected time limits', ...
    'Position',[.70 .55 .12 .5], ...
    'BackGroundColor',guiColor, ...
    'CallBack',@setLooping);
songThumbnailH=line( ...
    'Parent',sliderAxes, ...
    'XData',[0,0], ...
    'YData',[0,0], ...
    'Color',[1,1,1]);
smallPointerH=line( ...
    'Parent',sliderAxes, ...
    'XData',[0,0], ...
    'YData',[0,1], ...
    'LineWidth',3, ...
    'Color',[1, .5, .5]);
z=zoom(aH);
setAxesZoomMotion(z,aH,'horizontal');

% distAxes=axes(...
%     'Parent',DistPanel, ...
%     'Units','normalized',...
%     'Position', [.2,0.2,0.7,.8], ...
%     'Xlim',[0,1], ...
%     'Ylim',[0,1], ...
%     'Xtick',[], ...
%     'Ytick',[], ...
%     'LineWidth', .00001, ...
%     'Color', guiColor);
% featureDistPatch=patch(...
%     'Parent',distAxes, ...
%     'YData',zeros(1,nBins+2), ...
%     'XData',[0:1/(nBins-1):1,1,0], ...
%     'FaceAlpha',.2, ...
%     'EdgeAlpha',.4, ...
%     'FaceColor','r', ...
%     'EdgeColor','r');
% songDistPatch=patch(...
%     'Parent',distAxes, ...
%     'YData',[0,0], ...
%     'XData',[0:1/(nBins-1):1,1,0], ...
%     'FaceAlpha',.2, ...
%     'EdgeAlpha',.4, ...
%     'FaceColor','g', ...
%     'EdgeColor','g');
% 
% text('Parent',distAxes, ...
%     'String','song', ...
%     'FontSize',8, ...
%     'Units','normalized', ...
%     'Position', [0,-.2], ...
%     'Color', [0,.5,0]);
% text('Parent',distAxes, ...
%     'String','all songs', ...
%     'FontSize',8, ...
%     'Units','normalized', ...
%     'Position', [.2,-.2], ...
%     'Color', [.5,0,0]);
% noDataText=text('Parent',distAxes, ...
%     'String','NO SONG DATA', ...
%     'Units','normalized', ...
%     'Position',[.5,.5], ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment','Middle', ...
%     'FontSize',12, ...
%     'FontWeight', 'bold', ...
%     'Color',[0,.5,0], ...
%     'Visible', 'off');



%vertRect=[.1,.15;.9,.15;.9,.85;.1,.85];

%vrCreate=[0;1;1;0];

for featureInd=1:nFeatures
    selectFeatureButton{featureInd}  =   uicontrol(...
        'Parent', FeaturePanel, ...
        'Style','CheckBox', ...
        'Units','normalized',...
        'HandleVisibility','callback', ...
        'TooltipString',[sprintf('%s/',feat_info.fields{featureInd}{1:(end-1)}),regexprep(feat_info.names{featureInd},'.*/','')], ...
        'CallBack', @selectFeatureCallback, ...
        'String',feat_info.names{featureInd}, ...
        'Tag',num2str(featureInd), ...
        'BackGroundColor',guiColor, ...
        'Position',[0 (nFeatures-featureInd)/nFeatures 1 1]);%, ...
    if feat_info.isSongLevel(featureInd), set(selectFeatureButton{featureInd},'ForegroundColor',[.5,.5,.5], ...
            'TooltipString',[get(selectFeatureButton{featureInd},'TooltipString'), ' (song-level feature, only distribution shown)']); end %,'Enable','off'); end
    
    extent=get(selectFeatureButton{featureInd},'Extent');
    set(selectFeatureButton{featureInd},'Position', [0,.94*((nFeatures+1)-featureInd)/nFeatures,1,1.3*extent(4)]);
end

selectSong();
% uistack(fig,'top');

if exist('logo','var')
    
    logoaxes=axes(...
        'Parent',MainPanel, ...
        'Units','normalized',...
        'Position', [.3,0,.4,1], ...
        'Xlim',[0,1], ...
        'Ylim',[0,1], ...
        'Xtick',[], ...
        'Ytick',[], ...
        'LineWidth', .1, ...
        'Visible','off', ...
        'Color', guiColor);
    %,'Border','tight'
    image(logo,'Parent',logoaxes);
    
    set(logoaxes,'visible','off');
end


%%

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
            sliderPosition=get(playheads{1},'XData');
            sliderPosition=mean(sliderPosition(1:2))/(xlim(2)-xlim(1));
        else
            sliderPosition=mean(sliderHLim);
        end
        
        sliderHLim([1,4])=max(0,min(sliderPosition-sliderWidth/2, 1-sliderWidth));
        sliderHLim([2,3])=min(1,sliderHLim(1)+sliderWidth);
        set(sliderH,'XData',sliderHLim);
        set(aH,'Xlim',xlim(1)+sliderHLim(1:2)*(xlim(2)-xlim(1)));
        for i=1:length(featureAxes)
            set(featureAxes{i},'Xlim',xlim(1)+sliderHLim(1:2)*(xlim(2)-xlim(1)));
        end

        drawnow
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
            for i=1:length(featureAxes)
                set(featureAxes{i},'Xlim',xlim(1)+sliderLim*(xlim(2)-xlim(1))); % slider at the center of mouse movement
            end
        end
        drawnow
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
            for i=1:length(featureAxes)
                set(featureAxes{i},'Xlim',xlim(1)+sliderLim*(xlim(2)-xlim(1))); % slider at the center of mouse movement
            end
        end
        drawnow
    end

%PLAY
    function playPausePlayer(hObject, eventdata)
        
        if get(loopingButton,'Value')
            sliderHLim=get(sliderH,'XData');
            CurrentSelection=round(sliderHLim([1,2])*player.TotalSamples);
            CurrentSelection(1)=max(1,CurrentSelection(1));
            CurrentSelection(2)=min(player.TotalSamples,CurrentSelection(2));
            
            if CurrentSample<CurrentSelection(1)
                CurrentSample=max(CurrentSelection(1),CurrentSample);
            end
            
        end
        
        if not(ishandle(fig))
            return
        end
        
        state=get(PlayPauseButton,'Tag');
        
        CurrentFrame=max(1,sum(frameSummary<(xlim(1)+CurrentSample/Fs),2));
        
        for featureInd=1:length(playheads)
            set(playheads{featureInd},'XData',getxdata(framePos{featureInd}(:,CurrentFrame(featureInd))));
        end
        
        set(smallPointerH,'XData',(CurrentSample-1)/player.TotalSamples*[1,1]);
        
        if strcmp(state,'play') && not(isplaying(player))
            set(PlayPauseButton,'Tag','pause','cdata',pauseIcon);
            play(player,CurrentSample);
        elseif strcmp(state,'pause') && isplaying(player)
            set(PlayPauseButton,'Tag','play','cdata',playIcon);
            pause(player);
        else
            %do nothing
        end
        drawnow
    end

%STOP
    function stopPlayer(hObject, eventdata)
        
        if not(ishandle(fig))
            return
        end
        stop(player);
        set(PlayPauseButton,'Tag','play','cdata',playIcon);
        CurrentSample=get(player,'CurrentSample');
        CurrentFrame(:)=1;
        for featureInd=1:length(playheads)
            
            set(playheads{featureInd},'XData',xlim([1,1,1,1]));
        end
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
            
            if strcmp(get(PlayPauseButton,'Tag'),'pause')
                pause(player);
            end
            
            
            
            %CurrentFrame=length(find(framePos(1,:)<=CurrentPointAxes(1,1)));
            %pointerX=[framePos(1,CurrentFrame),framePos(2,CurrentFrame),framePos(2,CurrentFrame),framePos(1,CurrentFrame)];
            %pointerAlpha=pointerAlphaDefault;
            
            %CurrentFrame=max(1,sum(frameSummary<(xlim(1)+CurrentPointAxes(1,1)),2));
            for featureInd=1:length(playheads)
                set(playheads{featureInd},'XData',getxdata(framePos{featureInd}(:,CurrentFrame(featureInd))));
            end
            set(fig, 'WindowButtonMotionFcn', @draggingFcn)
            set(fig, 'WindowButtonUpFcn', @stopDragFcn)
        elseif all(CurrentPointSlider(1,1:2)<=1) && all(CurrentPointSlider(1,1:2)>=0) %mouse clicked inside sliderAxes
            %disable sliding during playback if followPointerButton
            %selected
            if strcmp(get(PlayPauseButton,'Tag'),'pause') && get(followPointerButton,'Value')
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
                    for featureInd=1:length(featureAxes)
                        set(featureAxes{featureInd},'Xlim',xlim(1)+sliderLim*(xlim(2)-xlim(1))); % slider at the center of mouse click
                    end
                    set(smallPointerH,'XData',(player.CurrentSample-1)/player.TotalSamples*[1,1]);
                end
                pointOnSlider=CurrentPointSlider(1,1)-sliderLim(1);
                if isequal(axesLim(:),xlim(:))
                    set(fig, 'WindowButtonMotionFcn', @draggingFcn);
                    set(fig, 'WindowButtonUpFcn', @stopDragFcn);
                else
                    set(fig, 'WindowButtonMotionFcn', @slideAxes);
                    set(fig, 'WindowButtonUpFcn', @stopSlide);
                end
            end
            
        else
            
            return
        end
        
        drawnow
        
    end

    function draggingFcn(varargin)
        %used while mouse button is held down
        CurrentPoint= get(aH, 'CurrentPoint');
        
        if CurrentPoint(1,1)<xlim(1) || CurrentPoint(1,1)>xlim(2)% || CurrentPoint(1,2)<ylim(1) || CurrentPoint(1,2)>ylim(2)
            return
        end
        CurrentFrame=max(1,sum(frameSummary<(xlim(1)+CurrentPoint(1,1)),2));
        for featureInd=1:length(selectedFeatures)
            set(playheads{featureInd},'XData',getxdata(framePos{featureInd}(:,CurrentFrame(featureInd))));
        end
        set(smallPointerH,'XData',CurrentPoint(:,1));
        drawnow
    end

    function stopDragFcn(varargin)
        set(fig,'WindowButtonMotionFcn', '');
        CurrentPoint= get(aH, 'CurrentPoint');
        if CurrentPoint(1,1)<xlim(1) || CurrentPoint(1,1)>xlim(2)% || CurrentPoint(1,2)<ylim(1) || CurrentPoint(1,2)>ylim(2)
            return
        end
        
        if strcmp(get(PlayPauseButton,'Tag'),'pause')
            pause(player);
            
        end
        %get current sample from current frame
        CurrentSample=floor(player.TotalSamples * (CurrentPoint(1,1)-xlim(1))/(xlim(2)-xlim(1)));
        
        CurrentFrame=max(1,sum(frameSummary<(xlim(1)+CurrentPoint(1,1)),2));
        for featureInd=1:length(playheads)
            set(playheads{featureInd},'XData',getxdata(framePos{featureInd}(:,CurrentFrame(featureInd))));
            
            %CurrentSample=min(CurrentSample,max(1,floor(Fs*(framePos{featureInd}(1,CurrentFrame(featureInd))-xlim(1))));
        end
        set(smallPointerH,'XData',(CurrentPoint(1,1)-xlim(1))/(xlim(2)-xlim(1))*[1,1]);
        
        
        
        %%if play button is activated
        if strcmp(get(PlayPauseButton,'Tag'),'pause')
            play(player,CurrentSample);
            
        end
        drawnow
    end

    function setLooping(hObject, eventData)
        if get(hObject,'Value')
            set(followPointerButton,'Value',false);
        end
        
    end

    function followPointer(hObject, eventData)
        if get(hObject,'Value')
            set(loopingButton,'Value',false);
        end
        
    end

    function selectSong(varargin)
        % select a song from miraudio struct and plot it in aH. Update also
        % possible curve and peak data.
        
        if not(ishandle(fig))
            return
        end
        
        %prevent distracting pushes of buttons
        for fi=selectFeatureButton(feat_info.isSongLevel==0)
            set(fi{1},'Enable','off');
        end
        
        if isequal(class(player),'audioplayer')
            stopPlayer();  % No selection
            player=[];
            java.lang.Runtime.getRuntime.gc; %Java garbage collection
            
        else
            CurrentSample=1;
        end
        
        songInd=get(audioPopupmenuH, 'Value');
        song = miraudio(feat_info.songNames{songInd});
        
        %start and end in seconds
        xlim=get(song,'Pos');
        xlim = xlim{1}{1}([1 end]);
        ylim=get(aH,'YLim');
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
        player.StopFcn = 'stopPlaying';
        set(aH,'Xlim',xlim);
        for featureInd=1:length(featureAxes)
            set(featureAxes{featureInd},'Xlim',xlim);
        end
        set(sliderAxes,'Xlim',[0,1]);
        
        %songImage=mirgetdata(miraudio(song,'Sampling',downSampleRate));
        songImage=mirgetdata(song);
        %songImage=.5+.5*(songImage-mean(songImage))/max(abs(songImage));
        
        songImage=(songImage-min(songImage))/(max(songImage)-min(songImage)); %normalize to [0,1], mean=0
        
        %songImagePos=(0:length(songImage)-1)./downSampleRate;
        
        %could handle stereo wave -> take mean across channels?
        
        set(songThumbnailH,'XData',(0:length(songImage)-1)/(length(songImage)-1),'YData', songImage, 'Color', songColor);
        %set(songH,'XData',songImagePos,'YData',.75*(ylim(2)-ylim(1))*songImage+mean(ylim),'Color',songColor);
        %selectedFeatureInds=find(selectedFeatures~=2);  %2 means hidden feature
        CurrentSample=1;
        for featureInd=1:length(selectedFeatures)
            cla(featureAxes{featureInd});
            set(featureAxes{featureInd},'Visible','off');
        end
        featureAxes={};
        framePos={};
        playheads={};
        frameSummary=[];
        
        
        for selectedFeature=selectedFeatures
            selectFeature(songInd,selectedFeature);
        end
        %if length(selectedFeatures)>0
        %    showFeatureStats(selectedFeatures(end));
        %end
        
        %prevent distracting pushes of buttons
        for fi=selectFeatureButton(feat_info.isSongLevel==0)
            set(fi{1},'Enable','on');
        end
        
        drawnow
    end

    function selectFeatureCallback(hObject, eventData)
        %prevent distracting pushes of buttons
        for fi=selectFeatureButton
            set(fi{1},'Enable','off');
        end
        songInd=get(audioPopupmenuH, 'Value');
        featureState=get(hObject,'Value');
        selectedFeature=str2double(get(hObject,'Tag'));
        if feat_info.isSongLevel(selectedFeature), set(hObject,'Value',false); end
        if featureState
            if feat_info.isSongLevel(selectedFeature)==0,
                selectedFeatures=[selectedFeatures,repmat(selectedFeature,1,feat_info.numAxes(selectedFeature))];
                selectFeature(songInd, selectedFeature);
            end
            %showFeatureStats(selectedFeature);
        else
            %remove feature
            if feat_info.isSongLevel(selectedFeature)==0
                removed=find(selectedFeatures==selectedFeature);
                for r = removed
                    cla(featureAxes{r});
                    set(featureAxes{r},'Visible','off');
                end
                for i = 1:(length(featureAxes))
                    set(get(featureAxes{i},'XLabel'),'visible','off');                    
                end
                %set(playheads{removed},'Visible','off');
                featureAxes(removed)=[];
                framePos(removed)=[];
                playheads(removed)=[];
                CurrentFrame(removed)=[];
                selectedFeatures(removed)=[];
                frameSummary(removed,:)=[];
                frameSummary=frameSummary(:,any(isfinite(frameSummary),1));
                
                scaleAxes(length(selectedFeatures));
                if ~isempty(selectedFeatures)
                    set(get(featureAxes{end},'XLabel'),'visible','on');
                end                
            end
        end
        
        %prevent distracting pushes of buttons
        for fi=selectFeatureButton
            set(fi{1},'Enable','on');
        end
        drawnow
    end

    function selectFeature(songInd, selectedFeature)
        
        %TODO
        %incorporate everything within this if statement into a more
        %general approach of handling multi-axis features.
        if feat_info.numAxes(selectedFeature)>1 && feat_info.cellinds(selectedFeature)==0
            return
            
            tmp=['FEATURES',sprintf('.%s',feat_info.fields{selectedFeature}{1:end})];
            if strcmp(tmp(end),'.'), tmp=tmp(1:(end-1)); end
            displayer(eval(tmp),featureAxes{nAxes},songInd);            
            disp_figure = gcf;
            set(gcf,'visible','off');
            if ~isa(eval(tmp),'mirtemporal')
                    framePos_tmp=get(eval(tmp),'FramePos');
             else
                framePos_tmp=get(eval(tmp),'Time');
                framePos_tmp{songInd}{1} = [framePos_tmp{songInd}{1}';framePos_tmp{songInd}{1}'];
            end
            axis_children = get(disp_figure ,'Children');
            print_title=true;
            for i = 1:size(axis_children)
                axis_child = axis_children(i);
                ch = get(axis_child,'Children');
                %this might change in the future
                if isa(ch,'matlab.graphics.primitive.Text')
                   feat_sub_number = get(ch,'String');
                else
                    featureAxes{length(featureAxes)+1}=copy(handle(aH));
                    nAxes=length(featureAxes);
                    set(featureAxes{nAxes},'Parent',get(aH,'Parent'));
                    scaleAxes(nAxes);
                    axes(featureAxes{nAxes});
                    set(ch,'Parent',featureAxes{nAxes})
                    %set(featureAxes{nAxes},'Children',get(axis_child,'Children'))
                    
                    set(featureAxes{nAxes},'Visible','on','Xlim',get(aH,'Xlim'));%,'ButtonDownFcn', @startDragFcn);
                    if print_title
                        set(get(featureAxes{nAxes},'Title'),'String',feat_info.titles{selectedFeature});                                                            
                        print_title=false;
                    end
                    set(get(featureAxes{nAxes},'YLabel'),'String',feat_sub_number);
                    
                    if length(framePos_tmp{songInd})>1 %|| isequal(features.types{selectedFeature},'miraudio')
                        framePos(nAxes)={[NaN;NaN]};
                    elseif iscell(framePos_tmp{songInd})
                        framePos(nAxes)=framePos_tmp{songInd};
                    else
                        framePos(nAxes)=framePos_tmp(songInd);
                    end
                    
                    %set(get(featureAxes{nAxes},'Children'),'AlphaData',.5)
                    ttl=get(featureAxes{nAxes},'Title');
                    set(ttl,'String',regexprep(get(ttl,'String'),',.*',''));
                    xlab = get(featureAxes{nAxes},'xlabel');
                    set(xlab,'String',get(xlab,'String'));
                    playheads{nAxes}=copy(handle(pointerH));
                    ydata=get(featureAxes{nAxes},'ylim');
                    set(playheads{nAxes},'Parent',featureAxes{nAxes},'YData',getydata(ydata));
                    set(featureAxes{nAxes},'ylim',ydata);
                    
                    frameSummary(nAxes,1:size(framePos{nAxes},2))=mean(framePos{nAxes});
                    frameSummary(frameSummary==0)=inf;
                    CurrentFrame=max(1,sum(frameSummary<(xlim(1)+CurrentSample/Fs),2));
                    %if ~strcmp(feat_info.types{selectedFeature},'miraudio') && ~strcmp(feat_info.types{selectedFeature},'mirenvelope')
                    %    set(playheads{nAxes},'XData',getxdata(framePos{nAxes}(:,CurrentFrame(nAxes))),'Visible','on');
                    %else
                    %    set(playheads{nAxes},'Visible','off');
                    %end
                    drawnow                    
                end
            end
            %%if play button is activated
            if strcmp(get(PlayPauseButton,'Tag'),'pause')            
                play(player,CurrentSample);
            end
            close(disp_figure);
            return
        end
        
        if not(ishandle(fig))
            return
        end
        
        %%if play button is activated
        if strcmp(get(PlayPauseButton,'Tag'),'pause')
            pause(player);
        end
        
        %create new axis and place the feature
        featureAxes{length(featureAxes)+1}=copy(handle(aH));
        nAxes=length(featureAxes);
        set(featureAxes{nAxes},'Parent',get(aH,'Parent'));
        scaleAxes(nAxes);
        axes(featureAxes{nAxes});
        if feat_info.cellinds(selectedFeature)>0
            displayer(eval(['FEATURES',sprintf('.%s',feat_info.fields{selectedFeature}{1:end}),'{',num2str(feat_info.cellinds(selectedFeature)),'}']),featureAxes{nAxes},songInd);
            framePos_tmp=get(eval(['FEATURES',sprintf('.%s',feat_info.fields{selectedFeature}{1:end}),'{',num2str(feat_info.cellinds(selectedFeature)),'}']),'FramePos');
        elseif feat_info.cellinds(selectedFeature)==0
            tmp=['FEATURES',sprintf('.%s',feat_info.fields{selectedFeature}{1:end})];
            if strcmp(tmp(end),'.'), tmp=tmp(1:(end-1)); end
            
            
            if feat_info.isMirdata(selectedFeature)
                displayer(eval(tmp),featureAxes{nAxes},songInd);
                %framePos_tmp=get(eval(tmp),'FramePos');
                if ~isa(eval(tmp),'mirtemporal')
                    framePos_tmp=get(eval(tmp),'FramePos');
                else
                    framePos_tmp=get(eval(tmp),'Time');
                    framePos_tmp{songInd}{1} = [framePos_tmp{songInd}{1}';framePos_tmp{songInd}{1}'];
                end
            else                
                %displayer(eval(tmp),featureAxes{nAxes},songInd);
                framePos_tmp=eval([tmp,'.framepos']);
                plot(mean(framePos_tmp{songInd}),eval([tmp,'.data{songInd}']));
            end
        else
            error('Check the feature importing.');
        end
        if length(framePos_tmp{songInd})>1 %|| isequal(features.types{selectedFeature},'miraudio')
            framePos(nAxes)={[NaN;NaN]};
        elseif iscell(framePos_tmp{songInd})
            framePos(nAxes)=framePos_tmp{songInd};
        else
            framePos(nAxes)=framePos_tmp(songInd);
        end
        
        set(featureAxes{nAxes},'Visible','on','Xlim',get(aH,'Xlim'));%,'ButtonDownFcn', @startDragFcn);
        %set(get(featureAxes{nAxes},'Children'),'AlphaData',.5)
        ttl=get(featureAxes{nAxes},'Title');
        set(ttl,'String',regexprep(get(ttl,'String'),',.*',''));
        xlab = get(featureAxes{nAxes},'xlabel');
        set(xlab,'String',get(xlab,'String'));
        playheads{nAxes}=copy(handle(pointerH));
        ydata=get(featureAxes{nAxes},'ylim');
        set(playheads{nAxes},'Parent',featureAxes{nAxes},'YData',getydata(ydata));
        set(featureAxes{nAxes},'ylim',ydata);
        frameSummary(nAxes,1:size(framePos{nAxes},2))=mean(framePos{nAxes});
        frameSummary(frameSummary==0)=inf;
        CurrentFrame=max(1,sum(frameSummary<(xlim(1)+CurrentSample/Fs),2));
        
        %if ~strcmp(feat_info.types{selectedFeature},'miraudio') && ~strcmp(feat_info.types{selectedFeature},'mirenvelope')
        %    set(playheads{nAxes},'XData',getxdata(framePos{nAxes}(:,CurrentFrame(nAxes))),'Visible','on');
        %else
        %    set(playheads{nAxes},'Visible','off');
        %end
        
        drawnow
        %%if play button is activated
        if strcmp(get(PlayPauseButton,'Tag'),'pause')            
            play(player,CurrentSample);
        end
    end


    function res = getydata(ydata)
        res=[ydata([1,1,2,2]),ydata(2)-.1*(ydata(2)-ydata(1)),ydata([2,2,1])];
    end

    function res = getxdata(xdata)
        s=xdata(2)-xdata(1);
        m=mean(xdata);
        res=[xdata([1,2,2])',m+.05*s,m,m-.05*s,xdata([1,1])'];
    end

    function scaleAxes(nAxes)
        apos=get(aH,'Position');
        scaling=(1.5./(1+nAxes))^.25;
        scaling2 = 1.005-0.005*nAxes;
        for axisInd=1:nAxes;
            
            set(featureAxes{axisInd},'Position',[apos(1), apos(2)+scaling2*(nAxes-axisInd)/(nAxes),apos(3),scaling2*scaling*apos(4)*1/nAxes]);
            if axisInd<nAxes,
                xl=get(featureAxes{axisInd},'XLabel');
                set(xl,'visible','off');
                ydata=get(featureAxes{axisInd},'ylim');
                set(playheads{axisInd},'YData',getydata(ydata));
            else
                xl=get(featureAxes{axisInd},'XLabel');
                
            end
        end
        
        
    end


%     function showFeatureStats(featureInd_stats)
%         songInd=get(audioPopupmenuH, 'Value');
%         
%         
%         ticklabels{1}=num2str(feat_info.valueRange(featureInd_stats,1));
%         ticklabels{2}=num2str(mean(feat_info.valueRange(featureInd_stats,:)));
%         ticklabels{3}=num2str(feat_info.valueRange(featureInd_stats,2));
%         
%         if length(ticklabels{1})>4
%             ticklabels{1}=num2str(feat_info.valueRange(featureInd_stats,1),'%1.2e');
%         end
%         if length(ticklabels{2})>4
%             ticklabels{2}=num2str(mean(feat_info.valueRange(featureInd_stats,:)),'%1.2e');
%         end
%         if length(ticklabels{3})>4
%             ticklabels{3}=num2str(feat_info.valueRange(featureInd_stats,2),'%1.2e');
%         end
%         
%         if feat_info.emptysong(songInd)
%             set(noDataText,'Visible','on')
%         else
%             set(noDataText,'Visible','off')
%         end
%         set(DistPanel,'Title',upper(regexprep(feat_info.names{featureInd_stats},'.*/','')));
%         set(distAxes,'Xtick',[0,.5,1],'XTickLabel',ticklabels);
%         set(featureDistPatch,'YData',[feat_info.distribution(featureInd_stats,:),0,0]);
%         if ~isempty(feat_info.songDistributions{featureInd_stats})
%             set(songDistPatch,'YData',[feat_info.songDistributions{featureInd_stats}(songInd,:),0,0]/2);
%         end
%         drawnow
%     end

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
