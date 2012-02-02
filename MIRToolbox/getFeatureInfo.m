function features=getFeatureInfo(f,nBins,smoothingFactor, songs)

%features=struct('names',{}, ...
%    'fields',{}, ...
%    'songNames',{}, ...
%    'data',{}, ...
%    'valueRange',[], ...
%    'distribution',[], ...
%    'songDistributions',[], ...
%    'frameLength',[], ...
%    'framediff',[], ...
%    'frameStart',[], ...
%    'hasPeaks',[], ...
%    'isSongLevel',[]);
features={};

fieldBranch={};
fieldInd=0;
nFeatures=0;
maxFieldNameLength=3;
binEdges=0:1/(nBins-1):1; %when visualizing the feature distributions, map feature values in the range of 0 and 1
songNames={};

if isstruct(f)
    recursiveCheck(fieldBranch,f);
elseif isequal(class(f),'mirscalar')
    fieldInd=1;
    fieldName=get(f,'Title');
    getData(f);
elseif iscell(f) && isequal(class(f{1}),'mirscalar')
        fieldInd=1;
        fieldName=get(f{1},'Title');
        getData(f{1});
else
    error('Feature set must be struct or mirscalar variable.');
end

    function recursiveCheck(fieldBranch,f)
        
        if isstruct(f)
            fieldInd=fieldInd+1;
            
            newFieldNames = fieldnames(f);
            numberOfFields = length(newFieldNames);
        end
        
        for current = 1:numberOfFields
            % Field information
            fieldName = newFieldNames{current};
            newField = f.(fieldName);
            
            if iscell(newField)
                if isequal(class(newField{1}),'mirscalar') %reached the end of branch and found a mirscalar feature
                    getData(newField{1});
                    
                end
            elseif isequal(class(newField),'mirscalar') %reached the end of branch and found a mirscalar feature
                    getData(newField);
                    
                   
            elseif isstruct(newField)
                fieldBranch{fieldInd}=fieldName;
                recursiveCheck(fieldBranch,newField);
            else
                continue
            end
            
            
        end
        fieldInd=fieldInd-1;
    end
if isempty(features)
    error('No mirscalar type of feature vectors in the set.');
else
    features.songNames=songNames{1};
end




    function getData(scalarFeature)
        
        featureData_tmp=get(scalarFeature,'Data');
        if ~isempty(songs)
            featureData_tmp=featureData_tmp(songs);
        end
        
        
        %read the branch of fieldnames to get summarized
        fieldBranch{fieldInd}=fieldName;
        featureName='';
        for i=1:fieldInd-1
            featureName=strcat(featureName,fieldBranch{i}(1:min(end,maxFieldNameLength)),'/');
        end
        featureName=strcat(featureName,fieldBranch{fieldInd});
        
        nFeatures=nFeatures+1; %add feature
        
        if length(featureData_tmp{1}{1})==1
            features.isSongLevel(nFeatures)=1;
        else
            features.isSongLevel(nFeatures)=0;
        end
        
        songNames{nFeatures}=get(scalarFeature,'Name'); %just to check that each feature is extracted from the same song. (what if the first song feature is empty due to the feature characteristics?)
        if ~isempty(songs)
            songNames{nFeatures}=songNames{nFeatures}(songs);
        end
        if nFeatures>1 && ~isequal(songNames{1},songNames{nFeatures})
            error('%s: all features must relate to the same set of songs.',featureName);
        end
        
        
        features.fields{nFeatures}=fieldBranch;
        features.names{nFeatures}=featureName;
        
        %peaks...
        peakPos=get(scalarFeature,'PeakPos');
        peakPos=peakPos{1};
        peakPos=peakPos{1};
        
        %if feature hasn't got peaks detected
        if isempty(peakPos) %not empty
            features.hasPeaks(nFeatures)=0;
        else
            features.hasPeaks(nFeatures)=1;
        end
        
        %get value range of the feature (summarize its overal distribution)
        features.valueRange(nFeatures,1:2)=[Inf,-Inf];
        features.distribution(nFeatures,1:nBins)=zeros(1,nBins);
        features.emptysong = zeros(1,length(featureData_tmp));
        
        for song=1:length(featureData_tmp)
            tmp = featureData_tmp{song}{1};
            if iscell(featureData_tmp{song}{1})
                tmp2 = [];
                for i = 1:length(tmp)
                    tmp2 = [tmp2 tmp{i}];
                end
                tmp = tmp2;
            end
            if ~isempty(tmp) && ~all(isnan(tmp))
                features.minsong(nFeatures,song) = min(tmp);
                minValue = min(tmp);
                features.maxsong(nFeatures,song) = max(tmp);
                maxValue = max(tmp);
            else
                warning('%s, song %d: No feature extracted. Check if there was some error in feature extraction. Including an empty feature...',featureName, song);
                features.minsong(nFeatures,song) = NaN;
                features.maxsong(nFeatures,song) = NaN;
                features.emptysong(song) = 1;
                continue
            end
            
            features.valueRange(nFeatures,1:2)=[ min(features.minsong(nFeatures,song),features.valueRange(nFeatures,1)), max(features.maxsong(nFeatures,song),features.valueRange(nFeatures,2)) ];
            
            
            
            %only sufficiently static frame length and diff allowed
            resolution=10000;
            allowedDiscrepancy=.5; %seconds
            
            %features.framediff(nFeatures,song)=mode(round((framePos(1,2:end)-framePos(1,1:end-1))*resolution)/resolution);
            %if max(abs((framePos(1,:)-features.framediff(nFeatures,song)*(0:length(framePos)-1))))>allowedDiscrepancy
            %    warning('%s, song %d: Frames are not equidistant. Check if there was some error in feature extraction. Including still...',featureName, song);
            %end
            %
            % if any(framePos(2,:)-framePos(1,:)<=0)
            %    warning('%s, song %d: Bad frame lengths (zero or negative). Check if there was some error in feature extraction. Including still...',featureName, song);
            %end
            
            %features.frameLength(nFeatures,song)=mode(round((framePos(2,:)-framePos(1,:))*resolution)/resolution);
            
            %features.frameStart(nFeatures,song)=framePos(1,1);
            %features.frameEnd(nFeatures,song)=framePos(2,end);
            
            %if features.isSongLevel(nFeatures)
            %    features.framediff(nFeatures,song)=features.frameEnd(nFeatures,song)-features.frameStart(nFeatures,song);
            %    features.frameLength(nFeatures,song)=features.frameEnd(nFeatures,song)-features.frameStart(nFeatures,song);
            %end
            
            
            %frames with peaks
            if features.hasPeaks(nFeatures)
                peakFrames=get(scalarFeature,'PeakPos');
                features.peakPos{nFeatures}{song}=peakFrames{song}{1}{1};
                
                %features.peakPos{nFeatures}{song}=features.frameStart(nFeatures,song)+features.framediff(nFeatures,song)*(peakPos-1);
                
                %get peak strength relative to feature value span
                features.peaks{nFeatures}{song}=(featureData_tmp{song}{1}(features.peakPos{nFeatures}{song})-minValue)/(maxValue-minValue);
            end
            
        end
        
        if iscell(featureData_tmp{song}{1})
            features.songDistributions{nFeatures} = [];
        else
            for song=1:length(featureData_tmp)
                %compute distribution of in one song, related to
                %the featureValueRange
                features.songDistributions{nFeatures}(song,1:nBins)=medfilt1(histc((featureData_tmp{song}{1}-features.valueRange(nFeatures,1))/(features.valueRange(nFeatures,2)-features.valueRange(nFeatures,1)),binEdges),smoothingFactor); %histogram, values related to the featureValueRange

                if isequal(features.songDistributions{nFeatures}(song,1:nBins),zeros(1,nBins)) %if filtering was too harsh due to small number of feature values
                    features.songDistributions{nFeatures}(song,1:nBins)=histc((featureData_tmp{song}{1}-features.valueRange(nFeatures,1))/(features.valueRange(nFeatures,2)-features.valueRange(nFeatures,1)),binEdges);
                end

                features.distribution(nFeatures,1:nBins)=features.distribution(nFeatures,1:nBins)+features.songDistributions{nFeatures}(song,1:nBins);

            end
        
            %map feature distributions to [0,1]
            features.songDistributions{nFeatures}=features.songDistributions{nFeatures}./repmat(max(features.songDistributions{nFeatures},[],2),1,nBins);
            features.distribution(nFeatures,1:nBins)=features.distribution(nFeatures,1:nBins)./repmat(max(features.distribution(nFeatures,1:nBins),[],2),1,nBins);
        end
            
        features.data{nFeatures}=featureData_tmp;
        framePos=get(scalarFeature,'FramePos');
        if ~isempty(songs)
            framePos = framePos(songs);
        end
        features.framePos{nFeatures}=framePos;
    end

end