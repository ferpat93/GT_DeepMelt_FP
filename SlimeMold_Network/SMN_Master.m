
%%% -  SMN_Master: Code that goes thru all the folders - %%%

%% Read Ledger and set constants %%

A_int = 1:0.5:4; % Area Intervals for part 1 (relative to initial area)
ind_G = find(A_int==3);

DataFolder = 'D:\Slime_Mold_Network\Data\Fusion';
ResultsFolder = 'D:\Slime_Mold_Network\Results\Fusion';

LedgerName = 'Ledger.xlsx';
T = readtable(fullfile(DataFolder,LedgerName));

%DataFolder = '/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data';
%ResultsFolder = '/Users/lfp3/Dropbox (GaTech)/GT/Spring-20/SMN/Data';

Folders=GetSubfolders(DataFolder);

%% ~ PRE - FUSION ~ %%

for f=100:numel(Folders) % For each Folder of images (each may contain several dishes)
    
    disp(strcat('Started processing folder: ',Folders{f}))
    Outputfile = fullfile(DataFolder,Folders{f},'DgdgT.mat');
    if isfile(Outputfile)
        disp('Already computed!')
        continue
    end
    
    ROI_filename = 'ROI_file.csv';
    srcFiles_filename = 'srcFiles.mat';
    ImgType = '*.jpg';
    %Cases={'Control','Control','Glucose','Glucose','NaCl','NaCl'};
    TF=T(ismember(T.Folder,Folders{f}),:);
    
    % Check if ROI's have been found
    if ~isfile(fullfile(DataFolder,Folders{f},ROI_filename)) % If diesn't exist ask for it
        % Get ROI's
    end   
        
    ROI = csvread(fullfile(DataFolder,Folders{f},ROI_filename));      
    nE = size(ROI,1);
    
    if isfile(fullfile(DataFolder,Folders{f},srcFiles_filename)) % If exists load it, otherwise compute it.
        srcFiles = load(fullfile(DataFolder,Folders{f},srcFiles_filename));
        srcFiles = extractfield(srcFiles,'srcFiles');
        srcFiles = srcFiles{1};
    else 
        srcFiles = dir(fullfile(DataFolder,Folders{f},ImgType));
        [imageIndexes,t]=sortTime(srcFiles);
        srcFiles = addvars(struct2table(srcFiles),imageIndexes,(1:numel(srcFiles))','NewVariableNames',{'time_order','file_order'});
        srcFiles = sortrows(srcFiles, 'time_order');
        srcFiles = addvars(srcFiles,t,'NewVariableNames','t');
        
        save(fullfile(DataFolder,Folders{f},srcFiles_filename),'srcFiles','-v7.3');
    end
    
    t = srcFiles.t; % times of images
    dataDishes = cell(nE,1);
    DT = cell(nE,1); % time vs distance vectors.
    Graphs = cell(nE,2);
    
    for e=1:nE % Loop over experiments of the folder

        disp(strcat('Started processing experiment #',num2str(e)))            
        tic
        mask = createCircleMask(ROI(e,3));
        
        
        % Process First Image (t=0)
        path = fullfile(DataFolder,Folders{f},srcFiles.name{1});
        [BW,centroids]= getBinary(path,ROI(e,:),[]);
        nCells = size(BW,3);
        Ao=permute(sum(BW,[1 2]),[3,2,1]); % Initial Areas
        if nCells>1; Do = DistanceCells(BW); end
        
        %% Process Image at fusion (time known)
        
        ft = str2double(TF(ismember(TF.N_Petri_Dish,e),:).Time_fusion{:}); % Fusion Time from table
        ift = GetFusionTime(ft,t,fullfile(DataFolder,Folders{f}),srcFiles,centroids,ROI(e,:));
        
        path = fullfile(DataFolder,Folders{f},srcFiles.name{ift}); % Path of image at fusion
        [BW, ~] = getBinary(path,ROI(e,:),centroids);
        Afo=permute(sum(BW,[1 2]),[3,2,1]);
        
        clear AT
        AT=cell(nCells,1); % Ledger for datapoints
        ReachTimes = zeros(1,numel(A_int));
        dataCells = cell(nCells,1);
        
        baseTimes = round(linspace(1,ift,22));
        
        for c=1:nCells           
            AT{c} = [AT{c};0 1 Do;t(ift) Afo(c)/Ao(c) 0]; 
            
            %% Evenly Spaced Values
            for i=2:numel(baseTimes)-1
                path = fullfile(DataFolder,Folders{f},srcFiles.name{baseTimes(i)}); % Path of image at fusion
                [BW, ~] = getBinary(path,ROI(e,:),centroids);              
                
                if nCells==1 % Just One cell -> Save Area
                    AT{1}=[AT{1}; t(baseTimes(i)) sum(BW,'all')./Ao nan];
                elseif size(BW,3)==1  % Two Cells, Fused -> No Area and Zero distance
                    AT{c}=[AT{c}; t(baseTimes(i)) nan 0];
                else % Two Cells, Not fused -> Areas and Distance
                    Ax=permute(sum(BW,[1 2]),[3,2,1])./Ao;
                    D = DistanceCells(BW);
                    AT{c}=[AT{c}; t(baseTimes(i)) Ax(c) D];
                    if c==1 % save for seccond cell as well
                        AT{2}=[AT{2}; t(baseTimes(i)) Ax(2) D];
                    end   
                end
                    
            end
            
            %% Interpolated to Target Area
            for Ai = 2:numel(A_int)
                A = A_int(Ai); % Objective Area
                cA = 0;
                dt=t(ift); % Initial dist between bounds
                epsA = 0.075; % Admisible difference between computed and objective (%)
                epsT = 6; % Minutes - epsilon between frames
                
                while and(dt>epsT,(abs(A-cA)/A)>epsA) % missing dt calc
                    [it,dt,AT{c}] = find_it(A,AT{c},t);
                    
                    if it>numel(t) % if SM never reaches that point
                        break
                    elseif dt==0 % Already computed
                        cA = AT{c}(find(AT{c}(:,1)==t(it)),2);
                        continue
                    end
                    
                    path = fullfile(DataFolder,Folders{f},srcFiles.name{it});
                    [BW, ~] = getBinary(path,ROI(e,:),centroids);
                    
                    if nCells==1 % Just One cell -> Save Area
                        AT{1}=[AT{1};t(it) sum(BW,'all')./Ao nan];
                    elseif size(BW,3)==1  % Two Cells, Fused -> No Area and Zero distance
                        AT{c}=[AT{c}; t(it) nan 0];
                    else % Two Cells, Not fused -> Areas and Distance
                        Ax=permute(sum(BW,[1 2]),[3,2,1])./Ao;
                        D = DistanceCells(BW);
                        AT{c}=[AT{c}; t(it) Ax(c) D];
                        cA = Ax(c);
                        if c==1 % save for seccond cell as well
                            AT{2}=[AT{2}; t(it) Ax(2) D];
                        end   
                    end
                end             
                if ~isinf(it); ReachTimes(Ai) = t(it); end
            end  % End Find Reach Times   
            
            [dataCells{c},Graphs{e,c}] = BF_Area_Indexes(c,fullfile(DataFolder,Folders{f}),srcFiles,[A_int' ReachTimes'],ROI(e,:),Ao); %% Gather Data         
        end
        
        Data_e = addvars(dataCells{1},repmat(e,numel(A_int),1),'NewVariableNames','nDish');
        if nCells>1            
            Data_e = [Data_e; addvars(dataCells{2},repmat(e,numel(A_int),1),'NewVariableNames','nDish')];
        end
        
        dataDishes{e}=Data_e;
        if nCells>1
            DT{e} = unique([AT{1}(:,[1 3]); AT{2}(:,[1 3])],'rows');
        end
    end % For Experiments in folder
    
    Graphs=[repelem(Folders(f),nE)' num2cell(1:nE)' Graphs];
           
    Data_f = addvars(dataDishes{1},repelem(Folders{f},size(dataDishes{1},1),1),'NewVariableNames','folder');
    for e=1:nE          
        Data_f = [Data_f; addvars(dataDishes{e},repelem(Folders{f},size(dataDishes{e},1),1),'NewVariableNames','folder')];
    end
    
    save(fullfile(DataFolder,Folders{f},'G.mat'),'Graphs','-v7.3');
    save(fullfile(DataFolder,Folders{f},'PreFusion_indexes.mat'),'Data_f','-v7.3');
    save(fullfile(DataFolder,Folders{f},'DT.mat'),'DT','-v7.3');
end % For Folders

%save(fullfile(DataFolder,'Graphs.mat'),'GraphArray','-v7.3');
    
%% ~ POST - FUSION ~ %%

TW = 180; % Time window (before/after fusion)
dt = 15; % Time Interval 

disp(' -> Started Post-Fusion Analysis  <- ')

for f=1:numel(Folders) % For each Folder of images (each may contain several dishes)
    
    disp(strcat('Started processing folder: ',Folders{f}))
    Outputfile = fullfile(DataFolder,Folders{f},'AfterFusion_indexes_blabla.mat');
    if isfile(Outputfile)
        disp('Already computed!')
        continue
    end
    
    ROI_filename = 'ROI_file.csv';
    srcFiles_filename = 'srcFiles.mat';
    ImgType = '*.jpg';
    %Cases={'Control','Control','Glucose','Glucose','NaCl','NaCl'};
    TF=T(ismember(T.Folder,Folders{f}),:); %Smaller table with just experiments of folder 'f'
    
    % Check if ROI's have been found
    if ~(isfile(fullfile(DataFolder,Folders{f},ROI_filename))) % If doesn't exist ask for it
        % Get ROI's
    end   
        
    ROI = csvread(fullfile(DataFolder,Folders{f},ROI_filename));      
    nE = size(ROI,1);
    
    if (isfile(fullfile(DataFolder,Folders{f},srcFiles_filename))) % If exists load it, otherwise compute it.
        srcFiles = load(fullfile(DataFolder,Folders{f},srcFiles_filename));
        srcFiles = extractfield(srcFiles,'srcFiles');
        srcFiles = srcFiles{1};
    else 
        srcFiles = dir(fullfile(DataFolder,Folders{f},ImgType));
        [imageIndexes,t]=sortTime(srcFiles);
        srcFiles = addvars(struct2table(srcFiles),imageIndexes,(1:numel(srcFiles))','NewVariableNames',{'time_order','file_order'});
        srcFiles = sortrows(srcFiles, 'time_order');
        srcFiles = addvars(srcFiles,t,'NewVariableNames','t');
        
        save(fullfile(DataFolder,Folders{f},srcFiles_filename),'srcFiles','-v7.3');
    end
    
    t = srcFiles.t; % times of images
    dataDishes = cell(nE,1);
    
    FusionRegionsExperiments = cell(nE,4);
    
    for e=1:nE % Loop over experiments of the folder

        disp(strcat('Started processing experiment #',num2str(e)))            
        tic
        mask = createCircleMask(ROI(e,3));

        % Process First Image (t=0) -> Just to get centroids and nCells
        % ((COMMON TO BOTH ANALYSES)) %
        FolderPath = fullfile(DataFolder,Folders{f});
        path = fullfile(FolderPath,srcFiles.name{1});
        [BW,centroids]= getBinary(path,ROI(e,:),[]);
        nCells = size(BW,3);
        Ao=permute(sum(BW,[1 2]),[3,2,1]); % Initial Areas
        
        %% -FUSION REGION ANALYSIS- %%
        
        % Process Image at fusion (time known)        
        ft = str2double(TF(ismember(TF.N_Petri_Dish,e),:).Time_fusion{:}); % Fusion Time from table
        
        if or(nCells<2,isnan(ft)) % No fusion!
            FusionRegionsExperiments{e,1} = [];
            FusionRegionsExperiments{e,2} = [];
            FusionRegionsExperiments{e,3} = [];
            FusionRegionsExperiments{e,4} = [];
            continue
        end
        
        [~,idx]=min(abs(t-(ft:dt:ft+TW)),[],1); % Get indexes of images at intervals requested (tf:tf+2h)              
        Fi = idx(1); % Reported Fusion TIME
        
        [FI, ~] = getBinary(fullfile(DataFolder,Folders{f},srcFiles.name{Fi}),ROI(e,:),centroids); %BF
        nc = size(FI,3);
        
        if nc>1 % Trouble!
            disp('reported fusion time seems disconnected')       
        end
        
        %sign = (nc > 1)*(1) + (~(nc > 1))*(-1);
        sign=-1;
        while size(FI,3)==1%abs(nc-size(FI,3))<1
            Fi=Fi+sign*1;
            [FI, ~] = getBinary(fullfile(DataFolder,Folders{f},srcFiles.name{Fi}),ROI(e,:),centroids);
        end

        %Fi = Fi -(size(FI,3)<2); % Find first image when the algorithm identified fusion (before observed fusion)
        [BF, ~] = getBinary(fullfile(DataFolder,Folders{f},srcFiles.name{Fi}),ROI(e,:),centroids); % Right before fusion
        %se=strel('disk',2);
        %CH=imdilate(bwconvhull(BF(:,:,2))+bwconvhull(BF(:,:,1)),se); % convex hull - 2 separate cells
        BFT=sum(BF,3); % Sum of cells before fusion  
        
        %[FF, ~] = getBinary(fullfile(DataFolder,Folders{f},srcFiles.name{idx(1)}),ROI(e,:),centroids); % At reported fusion
        [FF, ~] = getBinary(fullfile(DataFolder,Folders{f},srcFiles.name{Fi+5}),ROI(e,:),centroids); % Just after algorithm fusion time
        
        Dist=(bwdist(BF(:,:,1))+bwdist(BF(:,:,2))).*((bwdist(BF(:,:,1)).*bwdist(BF(:,:,2)))>0); % Distance between cells along fusion image
        Dist(Dist==0)=inf;
        minDist = min(Dist,[],'all');
        
        Diff = (FF-BFT)>0;
        
        [LL,nl] = bwlabel(Diff);
        CandReg = regionprops('table',LL, 'area', 'Centroid');
        CandReg = addvars(CandReg,(1:nl)','NewVariableNames','idx');
        CandReg(CandReg.Area<(minDist-1)/2,:)=[];
        CandReg = sortrows(addvars(CandReg,Dist(sub2ind(size(Dist),round(CandReg.Centroid(:,1)),round(CandReg.Centroid(:,2)))),'NewVariableNames','Dist'),'Dist');
  
        pcc = CandReg.idx;
        
        for ip = 1:numel(pcc)
            tempCC = bwconncomp(BFT+(LL==pcc(ip)));
            if tempCC.NumObjects ==1
                FusionArea=(LL==pcc(ip)); % Assumes there is only one region!
                tempRP = regionprops('table',FusionArea,'Centroid');
                centroid =tempRP.Centroid; % Centroid of fusion location
                break
            end
        end
        
        %ShortPath = (Dist<=(min(Dist.*FF,[],'all')+1)); % Pixels along the shortest path between cells
        % Difference between BF and FF
        %FusionArea = bwlabel(Diff.*ShortPath); % BW image with fusion area
        
%         nCCs =  max(FusionArea,[],'all'); % Check for more than 1 Fusion region:
%         if nCCs>1
%             disp('trouble!')
%         else
%             centroid =regionprops('table',FusionArea,'centroid').Centroid; % Centroid of fusion location
%         end
        
        iFA = FusionArea;  % Very little pixels with just the contact
        
        [I, ~] = getBinary(fullfile(DataFolder,Folders{f},srcFiles.name{idx(1)}),ROI(e,:),centroids); % At reported fusion (manual)
        
        maxDist = (ROI(e,3)/90)*5; % circle of 10mm in diameter aprox
        FusionArea = bwareafilt((iFA+and((bwdist(iFA).*I)<maxDist,(bwdist(iFA).*I))>0),1);
        
        %% First option to get Fusion Area
%         se=strel('disk',4);  
%         Seg = imdilate(bwlabel(imerode(I,se)),se); % Segmented Image by regions
%         idL=unique(Seg.*FusionArea); % Regions at the fusion region
%         idL(idL==0)=[];
%         
%         if numel(idL)==0
%             ii=0;
%             while numel(idL)<2
%                 %ii=ii+1;
%                 %se2=strel('disk',2*ii); 
%                 FusionArea = imdilate(FusionArea,se).*I;
%                 idL=unique(Seg.*FusionArea); % Regions at the fusion region
%                 idL(idL==0)=[];
%             end
%         end
%         
%         Region=zeros(size(FusionArea));
%         for j=1:numel(idL)
%             Region = Region + (Seg==idL(j)); % All regions at fusion
%         end   
% 
%       FusionArea = Region.*I; % Initial fusion Area (using observed fusion time)
             
        se=strel('disk',15);
        RegionOfAnalysis = imdilate(FusionArea,se); % Buffer of Analysis
        
        tempRP = regionprops('table',RegionOfAnalysis,'BoundingBox');
        Box =tempRP.BoundingBox; % Bounding box to crop image
        pad=10;
        Bounds = [ max(ceil(Box(2))-pad,1)  min(ceil(Box(2))+pad+Box(4),size(FusionArea,1)) ...
                max(ceil(Box(1))-pad,1)  min(ceil(Box(1))+pad+Box(3),size(FusionArea,2))];
        RegionOfAnalysis = RegionOfAnalysis(Bounds(1):Bounds(2),Bounds(3):Bounds(4)); % TRIM REGION OF ANALYSIS
        
        FusionRegions = zeros([size(RegionOfAnalysis) numel(idx)]);  
        FusionRegions(:,:,1)=FusionArea(Bounds(1):Bounds(2),Bounds(3):Bounds(4));
        
        Graphs = cell(numel(idx),1);
        Graphs{1} = getGraph_AF(I,centroids);
        
        for k=2:numel(idx) % For the rest of the times -> Get fusion area
            [I, ~] = getBinary(fullfile(DataFolder,Folders{f},srcFiles.name{idx(k)}),ROI(e,:),centroids); % Binary image at time intervals
            FusionRegions(:,:,k)=RegionOfAnalysis.*I(Bounds(1):Bounds(2),Bounds(3):Bounds(4)); % Pixels inside buffer
            Graphs{k} = getGraph_AF(I,centroids);
        end
        
        %% Store FusionRegions
        FusionRegionsExperiments{e,1} = FusionRegions;
        FusionRegionsExperiments{e,2} = Bounds;
        FusionRegionsExperiments{e,3} = RegionOfAnalysis;  
        FusionRegionsExperiments{e,4} = Graphs;  

        %% -GENERAL INDEXES- %% -> Before And After fusion 
        % Single table with 5 points per experiment: -2h -1h 0 1h 2h -> Add
        Tgi = [(-1:0.5:1)'.*TW (linspace(ft-TW,ft+TW,5))']; % times general indexes
       
        % [Cell#, Area Normalized, time, initialArea, TotalArea, perc.
        % Clusters, Area Enclosed, Area WhiteSpace, #WS regions, Width dist,
        % Ldist]
    
        dataDishes{e}=AF_Area_Indexes(e,srcFiles,Tgi,ROI(e,:),Ao,centroids,FolderPath);     
        
    end % Loop Experiments of folder
      
    Data_f = addvars(dataDishes{1},repelem(Folders{f},size(dataDishes{1},1),1),'NewVariableNames','folder');
    for e=2:nE         
        if ~isempty(dataDishes{e})
            Data_f = [Data_f; addvars(dataDishes{e},repelem(Folders{f},size(dataDishes{e},1),1),'NewVariableNames','folder')];
        end
    end
    
    % Save General Indexes
    save(fullfile(DataFolder,Folders{f},'AfterFusion_indexes.mat'),'Data_f','-v7.3');      
    
    %Save Fusion Regions
    save(fullfile(DataFolder,Folders{f},'FusionRegions.mat'),'FusionRegionsExperiments','-v7.3');    
        
end % For Folders

%% AUXILIAR FUNCTIONS

function [names]=GetSubfolders(path)
        d = dir(path);
        isub = [d(:).isdir]; %# returns logical vector
        names = {d(isub).name}';
        names(ismember(names,{'.','..'})) = [];
end 
function [indexes,t]=sortTime(srcFiles)
    datetomin=6.94444e-4;
    DataFolder=srcFiles(1).folder;
    nFiles=numel(srcFiles);
    dates=cell(nFiles,1);
    for i=1:nFiles
      dates{i}=extractfield(imfinfo(fullfile(DataFolder,srcFiles(i).name)),'DateTime');  
    end
    
    d=cell2mat(cellfun(@(x) datenum(x,'yyyy:mm:dd HH:MM:SS')-2019*365,dates,'UniformOutput',false));
    [~,idx]=sort(d);
    DD=sortrows([idx (1:numel(d))'],1);
    indexes=DD(:,2);
    t=cumsum(round((d(idx(2:end))-d(idx(1:end-1)))./datetomin));
    t=[0;t];
end
function [mask] = createCircleMask(L)
    
    [xx,yy] = meshgrid(1:L,1:L);
    mask = false(L,L);
    mask = mask | hypot(xx - L/2, yy - L/2) <= 0.99*(L/2);

end
function [BW,centroids] = getBinary(path,ROI,centroids)

%% Load Image
    Io = imread(path);
    Io = Io(ROI(2):ROI(2)+ROI(3)-1,ROI(1):ROI(1)+ROI(3)-1,:);
    mask = createCircleMask(ROI(3));
    
%% Preprocess image -> Raw to Full Binary

    % pixels per mm
    D_dish = 90; %mm
    ppmm = size(Io,1)/D_dish;
    
    %%% Enhance contrast - grayscale - binarize

    % Enhance
    sigma = 0.4;
    alpha = 0.675;
    I = locallapfilt(Io, sigma, alpha); % Laplacian filtering -  

    lab = rgb2lab(I); % to LAB
    N=lab(:,:,1)-lab(:,:,3); % Compound 
    GI=(N-min(N,[],'all'))/(max(N,[],'all')-min(N,[],'all')); % Normalize
    GI = imgaussfilt(GI,2); % Apply gaussian filter to blur
    
    % Binarize
    BI = (~imbinarize(GI,'adaptive','ForegroundPolarity','dark','Sensitivity',0.55).*(mask));
    
    se = strel('disk',25);
    GI2 = imcomplement(imtophat(imcomplement(GI),se));

    wsi=3/100; % Watershed value - values around 3 and 6% -> depends on gaussian filter

    % Skeleton from watershed
    skelImage = (watershed(imhmin(imcomplement(GI2),wsi))==0);

    BW = ((BI+skelImage)>0).*mask; % Full binary Image  
    
%% Split Cells and measure indexes
    
    T = ppmm^2 * 50; %threshold value to discard noise: considers anything larger than 30 mm2
          
    if isempty(centroids) % Get centroids if its the 1st frame
        BL = bwlabel(BW);
        numPixels = regionprops('table',BL,'Area','Centroid');
        idx = find(numPixels.Area>T);  
        NC = numel(idx);       
        if (NC>2)
            warning('More than 2 large regions in first image')
            pauskjkjke
        else
            centroids=numPixels.Centroid(idx,:);
        end
    end 
    
    nCells = size(centroids,1); % 1 or 2 initial cells
    
    %% Split cells in two
    
    % Connected Print %
    %BW = ConnectedPrint(BW,T,nCells);   
    BW = ConnectedPrint2(BW,T,centroids);   
    
    CC = bwconncomp(BW);
    NC=CC.NumObjects;
    
    if NC<=nCells
        
        DC=inf(nCells,NC); %Rows are centroids, Cols are CC's
        BW=zeros(size(BW,1),size(BW,2),NC);
        
        for k=1:NC
            [row,col] = ind2sub(size(I),CC.PixelIdxList{k});
            for cs=1:nCells
                DC(cs,k)=min(min(((row-centroids(cs,2)).^2+(col-centroids(cs,1)).^2).^(1/2)),DC(cs,k));
            end
            c=zeros(size(BI));
            c(CC.PixelIdxList{k}) = 1;
            BW(:,:,k)=c;
        end
    
        [m,In]=min(DC,[],2);
        if max(m)>25
            warning(strcat('Ojo! Far away centroids - image: ',num2str(ind_im)))
        elseif and(NC==2,nCells==2)
            if In(2)-In(1)==0
                warning('Ojo! not one to one relation')   
            end
        end
        
    else
        
        warning(strcat('More CCs than initial cells!',num2str(ind_im))) 
        pause
    end
    
end
function [Connected2]=ConnectedPrint(BW,T,nCells)
    
    BL = bwlabel(BW);
    numPixels = regionprops('table',BL,'Area');
    regionsArea = sortrows([ numPixels.Area (1:numel(numPixels))'],'descend');
    NC = sum(regionsArea(:,1)>T); % Number of cells
    
    HNC = min(nCells,NC); % Hypotesis on the number of cells
    Connected = zeros(size(BW));
    for i=1:HNC
        Connected = Connected + (BL==regionsArea(i,2));
    end
    
    DistConnected=bwdist(Connected); 
    
    %RegionsToAdd=regionsArea(find(and(regionsArea(:,1)<T,regionsArea(:,1)>15)),:); % Leftover regions larger than 15px but less than T
    maxDist = 25; % Max dist to add  
   
    %% Add small Areas
    Connected2=Connected;
    for i=HNC+1:size(regionsArea,1)
        if regionsArea(i,1)>15
            D=DistConnected+bwdist(BL==regionsArea(i,2));  % Distance between Connected and region i
            if min(D,[],'all')/2<=maxDist % If its close enough
                Connected2 = (Connected2 + (BL==regionsArea(i,2)) + (D<=1+min(D,[],'all')))>0;
            end
        end
    end
    
    %% Check the new areas are connected to large ones
    
    [BL,NC2] = bwlabel(Connected2);

    if NC2>NC
        
        numPixels = regionprops('table',BL,'Area');
        regionsArea = sortrows([ numPixels.Area (1:numel(numPixels))'],'descend');
        NC2 = sum(regionsArea(:,1)>T); % Number of cells
        
        Connected2 = zeros(size(BW));
        for i=1:NC2
            Connected2 = Connected2 + (BL==regionsArea(i,2));
        end
    end
    
    % (at this point we know NC2<=NC)
    
    if NC2>nCells %  trouble

        DistConnected=bwdist(Connected2); 
        LL = bwlabel(Connected2);      
        %Connected2 = zeros(size(BW));
        
        for i=nCells+1:NC2
            D=DistConnected + bwdist(LL==i);  % Distance between Connected and region i
            Connected2 = (Connected2 + (LL==i) + (D<=2+min(D,[],'all')))>0;
        end  

    end   
        
end
function [it,dt,AT] = find_it(A,AT,times) % Fragile when not 1-to-1! -> Then deadlocks (not unique)
    % AT: Col 1 - Time /*/ Col 2 - Areas (normalized) 
    
    AT = unique(sortrows(AT,1),'rows'); % Sort by time
    low = find((AT(:,2)<A),1,'last');
    up = find((1:size(AT,1)>low)'.*(~isnan(AT(:,2))),1); % Find index of first image with larger area.
    
    if or(isempty(low),isempty(up)) % If SM never reaches the given area
        it=inf;
        dt=inf;      
        return
    end
    
    dt = AT(up,1)-AT(low,1);
    dA = AT(up,2)-AT(low,2);    
    m= dt/dA;    
    t=(A-AT(low,2))*m+AT(low,1);
 
    dist = abs(times - t);
    it = find(dist == min(dist));
    
    if (any(AT(:,1)==times(it)))% Repeated location
        dt=0;
    end
    
end
function[ift] = GetFusionTime(ft,Ts,folderPath,srcFiles,centroids,ROI)

    nCells = size(centroids,1);
    eps = 10 ; % Epsilon: 10 min of difference
    
    if or(nCells<2,isnan(ft)) % No Fusion at all! Find peak
        ft = round(Ts(end)*0.75); % Assume peak is at 3/4 of the total time              
        dist = abs(Ts - ft);
        ift = find(dist == min(dist)); % Index fusion time 
        
        % Leave like this for now, to add the peak finding section
          
    else % Two Cells and fusion is observed, find index
        
        dist = abs(Ts - ft);
        ift = find(dist == min(dist),1); % Index fusion time - minus 

        path = fullfile(folderPath,srcFiles.name{ift});
        [BW,~]= getBinary(path,ROI,centroids);
        nComp = size(BW,3);
        changed = false;
    
        if (nComp==1) % fusion already occurred try previous images
            sign = -1;
        else % fusion not yet ocurred -> try posterior images
            sign = 1;
        end  
        
        ii=ift;
        
        while ~changed % Do until nComp changes

            fi = round(ii + sign*(eps/(Ts(ii+1)-Ts(ii))));         
            path = fullfile(folderPath,srcFiles.name{fi});
            [BW,~]= getBinary(path,ROI,centroids);        
            
            if abs(nComp-size(BW,3))>0.1 % if changed
                changed = ~changed;             
            else
                ii = fi;
            end                    
        end
        
        if (nComp==1) % fusion already occurred try previous images
            ift = fi;
        else % fusion not yet ocurred -> try posterior images
            ift = ii;
        end  
        
        
    end


    
end
function [data, G] = BF_Area_Indexes(c,folderPath,srcFiles,RT,ROI,Ao) % returns a table

    centroids = [];
    t = srcFiles.t;    
    nA = size(RT,1);

    nBins=20;
    Wbins=linspace(0,2.5,nBins+1); % Width Bins 
    Lbins=linspace(0,6,nBins+1); % Length Bins
    
    % [Cell#, Area Normalized, time, initialArea, TotalArea, perc.
    % Clusters, Area Enclosed, Area WhiteSpace, #WS regions, Width dist,
    % Ldist, W fit par, L fit par,numNodes,numEdges, Length Paths, Length Euclidean, meanNodeDegree, alpha]
    TableVariables = {'cell','An','t','Ao','At','pC','Ae','Aws','nwsr','Wbins','Lbins','Wpar','Lpar',...
        'nNodes_nEdges','Le_p'};
    AreaInds = zeros(size(RT,1),6);
    W = zeros(nA,nBins);
    L = zeros(nA,nBins);
    Wpar = zeros(nA,2);
    Lpar = zeros(nA,2);
    inds = zeros(nA,4);
    
    for it = 1:nA % get indexes for each Area
        
        if  RT(it,2)>t(end)+1 % if SM never reaches that point           
            AreaInds(it,:) = nan(1,6);
            W(it,:) = nan(1,nBins);
            L(it,:) = nan(1,nBins);
            continue
        end
        
        dist = abs(t - RT(it,2));
        ift = find(dist == min(dist),1); % Index fusion time 
        
        path = fullfile(folderPath,srcFiles.name{ift});
        [BW,centroids]= getBinary(path,ROI,centroids); % Corresponding Binary Image 

        BW = BW(:,:,c); % Get only the one we need
        centroid = centroids(c,:);
        
        %% Split BW into Veins/Clusters
        
        % Identify Clusters and Measure Vein Width
        D_dish = 90; %mm
        ppmm = size(BW,1)/D_dish;   
        maxVeinSize = round(1*ppmm)+1; % 4mm meaning 8mm min dimension.
         
        seE = strel('disk', double(maxVeinSize)); % erode and dilate image
        EI = imerode(BW,seE);
        seD = strel('disk', (round(1.5*maxVeinSize)+1));
        DI = imdilate(EI,seD);

        clustersI = bwlabel(DI.*(BW)); % Get Clusters
        clusters_label = unique(clustersI.*EI);
        clustersI = imfill(ismember(clustersI, clusters_label(2:end)),'holes');

        Full_BI = (BW + clustersI)>0; % after filling inside of cluster
        Idist = bwdist(~Full_BI); % Binary distance (without clusters!)
        Veins = Idist.*(~clustersI);

        fullBW = logical((clustersI + Veins>0)>0);
        filledBW = imfill(fullBW,'holes');
        whitespace = filledBW-fullBW;
        
        total_Area = sum(fullBW,[1 2]); % Total SM Area
        p_clusters = 100*(sum(clustersI,[1 2])/total_Area); % Percentage of Area that is Clusters
        white_area = sum(whitespace,[1 2]); % Total Whitespace Area (trapped area)
        nWSA = max(bwlabel(whitespace),[],'all'); % Number of whistespace regions
        
        AreaInds(it,:) = [Ao(c) total_Area p_clusters total_Area+white_area white_area nWSA];
        
        %% Distribution of W and L
        
        skel_refined = bwskel(Full_BI>0,'MinBranchLength',10);
        skel_dist = skel_refined.*Idist;
        
        if it==nA; needG=1; else; needG=0; end
        [W(it,:),Wpar(it,:),L(it,:),Lpar(it,:),inds(it,:),G] = VeinsDistribution(skel_dist,clustersI,Wbins,Lbins,needG,centroid); %% Call to external function      
        
    end
    
    data = array2table([repmat(c,size(RT,1),1) RT AreaInds],'VariableNames',TableVariables(1:9));
    data = addvars(data,W,L,Wpar,Lpar,inds(:,1:2),inds(:,3:4),'NewVariableNames',TableVariables(10:end));
         
end
function [data] = AF_Area_Indexes(e,srcFiles,RT,ROI,Ao,centroids,folder) % returns a table

    t = srcFiles.t;    
    nA = size(RT,1);

    nBins=20;
    Wbins=linspace(0,2.5,nBins+1); % Width Bins 
    Lbins=linspace(0,6,nBins+1); % Length Bins
    
    % [Cell#, Area Normalized, time, initialArea, TotalArea, perc.
    % Clusters, Area Enclosed, Area WhiteSpace, #WS regions, Width dist,
    % Ldist]
    TableVariables = {'nDish','tfusion','t','Ao','At','pC','Ae','Aws','nwsr','Wbins','Lbins','Wpar','Lpar','Inds'};
    AreaInds = nan(nA,6,2);
    W = nan(nA,nBins,2);
    L = nan(nA,nBins,2);
    Lpar = nan(nA,2,2);
    Wpar = nan(nA,2,2);
    inds = nan(nA,4,2);
    
    for it = 1:nA % get indexes for each Area
        
        if  RT(it,2)>t(end)+1 % if SM never reaches that point           
            continue
        end
             
        dist = abs(t - RT(it,2));
        ift = find(dist == min(dist),1); % Index fusion time 
        
        path = fullfile(folder,srcFiles.name{ift});      
        [BWA,~]= getBinary(path,ROI,centroids); % Corresponding Binary Image 
        nCells = size(BWA,3); 
        
        for c=1:nCells
            
            BW = BWA(:,:,c); % Get only the one we need
            
            %% Split BW into Veins/Clusters

            % Identify Clusters and Measure Vein Width
            D_dish = 90; %mm
            ppmm = size(BW,1)/D_dish;   
            maxVeinSize = round(1*ppmm)+1; % 4mm meaning 8mm min dimension.

            seE = strel('disk', double(maxVeinSize)); % erode and dilate image
            EI = imerode(BW,seE);
            seD = strel('disk', (round(1.5*maxVeinSize)+1));
            DI = imdilate(EI,seD);

            clustersI = bwlabel(DI.*(BW)); % Get Clusters
            clusters_label = unique(clustersI.*EI);
            clustersI = imfill(ismember(clustersI, clusters_label(2:end)),'holes');
            labeled_clusters = int8(bwlabel(clustersI));

            Full_BI = (BW + clustersI)>0; % after filling inside of cluster
            Idist = bwdist(~Full_BI); % Binary distance (without clusters!)
            Veins = Idist.*(~clustersI);

            fullBW = logical((clustersI + Veins>0)>0);
            filledBW = imfill(fullBW,'holes');
            whitespace = filledBW-fullBW;

            total_Area = sum(fullBW,[1 2]); % Total SM Area
            p_clusters = 100*(sum(clustersI,[1 2])/total_Area); % Percentage of Area that is Clusters
            white_area = sum(whitespace,[1 2]); % Total Whitespace Area (trapped area)
            nWSA = max(bwlabel(whitespace),[],'all'); % Number of whistespace regions

            AreaInds(it,:,c) = [Ao(c) total_Area p_clusters total_Area+white_area white_area nWSA];

            %% Distribution of W and L

            skel_refined = bwskel(Full_BI>0,'MinBranchLength',10);
            skel_dist = skel_refined.*Idist; 
            [W(it,:,c),Wpar(it,:,c),L(it,:,c),Lpar(it,:,c),inds(it,:,c),~] = VeinsDistribution(skel_dist,labeled_clusters,Wbins,Lbins,0,[]); %% Call to external function  
        end
        
        if nCells ==2
            AreaInds = [ sum(AreaInds(:,1:2,:),3) sum(AreaInds(:,3,:).*AreaInds(:,2,:),3)./sum(AreaInds(:,2,:),3) sum(AreaInds(:,4:6,:),3)];
        else
            AreaInds = AreaInds(:,:,1);
        end
        
        W = nanmean(W,3);
        L = nanmean(L,3);
        Wpar = nanmean(Wpar,3);
        Lpar = nanmean(Lpar,3);
        inds = nanmean(inds,3);
        
    end
    
    data = array2table([repmat(e,size(RT,1),1) RT AreaInds],'VariableNames',TableVariables(1:9));
    data = addvars(data,W,L,Wpar,Lpar,inds,'NewVariableNames',TableVariables(end-4:end));
    
end
function [Connected2]=ConnectedPrint2(BW,T,centroids)
    
    nCells = size(centroids,1);
    centroids = round(centroids);
    
    BL = bwlabel(BW);
    numPixels = regionprops('table',BL,'Area');
    regionsArea = sortrows([ numPixels.Area (1:numel(numPixels))'],'descend');
    NCCs = sum(regionsArea(:,1)>T); % Number of cells
    
    DC=inf(nCells,NCCs); %Rows are centroids, Cols are CC's
        
    for cs=1:nCells 
        CI = zeros(size(BW));
        CI(centroids(cs,2),centroids(cs,1))=1;
        CID = bwdist(CI);
      
        for k=1:NCCs
            A = CID.*(BL==regionsArea(k,2));
            DC(cs,k)=min(A(A > 0),[],'all');
        end
    end
    
    [m,In]=min(DC,[],2);
    
    iCC = unique(In);
    NC = numel(iCC); % number of cells in the current step (<=nCells)
    
    if sum(m>25)>0
        disp(' Far Away from the centroids')
    end  

    Connected = zeros(size(BW));
    for i=1:NC 
        Connected = Connected + (BL==regionsArea(iCC(i),2));
    end
    
    % Add Left over regions
    DistConnected=bwdist(Connected);   
    %RegionsToAdd=regionsArea(find(and(regionsArea(:,1)<T,regionsArea(:,1)>15)),:); % Leftover regions larger than 15px but less than T
    maxDist = 25; % Max dist to add  
   
    %% Add small Areas
    Connected2=Connected;
    regionsArea(iCC,:)=[];
    for i=1:size(regionsArea,1)
        if regionsArea(i,1)>15
            D=DistConnected+bwdist(BL==regionsArea(i,2));  % Distance between Connected and region i
            if min(D,[],'all')/2<=maxDist % If its close enough
                Connected2 = (Connected2 + (BL==regionsArea(i,2)) + (D<=1+min(D,[],'all')))>0;
            end
        end
    end
    
    %% Check the new areas are connected to large ones
    
    [BL,NC2] = bwlabel(Connected2);

    if NC2>NC
        Connected2 = bwareafilt(Connected2,NC);
%         numPixels = regionprops('table',BL,'Area');
%         regionsArea = sortrows([ numPixels.Area (1:numel(numPixels))'],'descend');
%         NC2 = sum(regionsArea(:,1)>T); % Number of cells
%         
%         Connected2 = zeros(size(BW));
%         for i=1:NC2
%             Connected2 = Connected2 + (BL==regionsArea(i,2));
%         end
    end
    
    % (at this point we know NC2<=NC)
    
    if NC>nCells %  trouble -> shouldn't be here
        
        DistConnected=bwdist(Connected2); 
        LL = bwlabel(Connected2);      
        %Connected2 = zeros(size(BW));
        
        for i=nCells+1:NC2
            D=DistConnected + bwdist(LL==i);  % Distance between Connected and region i
            Connected2 = (Connected2 + (LL==i) + (D<=2+min(D,[],'all')))>0;
        end  

    end   
        
end
function [D] = DistanceCells(BF) % Only called if nCells ==2
        px1mm = size(BF,1)/90; % number of pixels to 1 mm
        Dist=(bwdist(BF(:,:,1))+bwdist(BF(:,:,2))).*((~BF(:,:,1).*~BF(:,:,2))>0); % Distance between cells along fusion image
        Dist(Dist==0)=inf;
        D = min(Dist,[],'all')/(2*px1mm);
end

function [G] = getGraph_AF(BW,centroids)

    %% Split BW into Veins/Clusters

    % Identify Clusters and Measure Vein Width
    D_dish = 90; %mm
    ppmm = size(BW,1)/D_dish;   
    maxVeinSize = round(1*ppmm)+1; % 4mm meaning 8mm min dimension.

    seE = strel('disk', double(maxVeinSize)); % erode and dilate image
    EI = imerode(BW,seE);
    seD = strel('disk', (round(1.5*maxVeinSize)+1));
    DI = imdilate(EI,seD);

    clustersI = bwlabel(DI.*(BW)); % Get Clusters
    clusters_label = unique(clustersI.*EI);
    clustersI = imfill(ismember(clustersI, clusters_label(2:end)),'holes');

    Full_BI = (BW + clustersI)>0; % after filling inside of cluster
    Idist = bwdist(~Full_BI); % Binary distance (without clusters!)
    skel_refined = bwskel(Full_BI>0,'MinBranchLength',10);
    skel_dist = skel_refined.*Idist;

    % Build Graph
    skel = skel_dist>0;
    px1mm = size(skel,1)/90; % number of pixels to 1 mm
    
    G=[]; % Initialize empty graph
    
    %% Check For empty graph
    if sum((skel-clustersI)>0,'all')<5  % Just the Cluster! - Create graph with a single node.
        disp('Empty graph!')     
        return
    end
       
    % I is the binary skeleton to morph into a graph
    BP = bwmorph(skel,'branchpoints');
    EP = bwmorph(skel,'endpoints');
    
    I_Nodes = EP + BP;
    nNodes = sum(I_Nodes,'all');
    
    % Get Dilated Node Image - 3
    SE = strel('square', 3);
    LN3 = imdilate(I_Nodes,SE);
 
    %% EDGES: ['Width' 'NumPixels']
    branches = (skel - (LN3>0))>0;
    LB = bwlabel(branches);  
    nEdges = max(LB(:));

    Per = regionprops(LB,'Perimeter');
    endpts = bwmorph(branches,'endpoints').*LB;
      
    if nEdges>0       
        Edges=zeros(nEdges,3); % Thickness, Le, Lp
        
        for b = 1:nEdges

            % Set Edge thickness
            Ds = skel_dist(LB==b);
            Ds(Ds==1)=[];
            if isnan(mean(Ds))
                D = 2.25;
            else
                D = 2*mean(Ds);
            end

            [rows,cols] = find(endpts==b);
            
            if numel(rows) ~= 2
              Le = NaN;
              Lp= nan;
            else
              Le = sqrt((rows(1)-rows(2))^2+(cols(1)-cols(2))^2);
              Lp = Per(b).Perimeter/2;
            end
          
            Edges(b,:) = [D Le Lp]./px1mm;

        end
     
        [LN,nNodes] = bwlabel(I_Nodes);
        SE = strel('square', 3);
        LN3 = imdilate(LN,SE);

        SED = strel('disk', 2);
        LB2 = imdilate(LB,SED);

        %% NODES : Graph Nodes: [ X Y Type* ]    *:comes later  
        labeled_clusters = int8(bwlabel(clustersI));
        s = regionprops(labeled_clusters,'centroid');
        [~,ind]=min(pdist2(centroids,cat(1,s.Centroid)),[],2);
        clustersI = clustersI + (labeled_clusters==ind(1))+ 2.*(labeled_clusters==ind(2));

        GraphNodes=zeros(nNodes,3);
        GraphNodes(:,1:2) = table2array(regionprops('table',LN,'centroid'));
        GraphNodes(:,3) = clustersI(sub2ind(size(LN),round(GraphNodes(:,2)),round(GraphNodes(:,1))));

        nEdges = size(Edges,1);
        GraphEdges = [zeros(nEdges,2) Edges];

        for b = 1:nEdges
            % Get Parent Nodes
            nb = unique(LN3(LB2==b))';
            nb(nb==0)=[];       
            if numel(nb)>2
                nb = nb(1:2); % Really bad way to do it -> OJO!
            elseif numel(nb)<2
                nb = [nb nb];
            end      
            GraphEdges(b,1:2)=nb;
        end   
    
        %% Set Graph  
        EdgeTable = table(GraphEdges(:,1:2),GraphEdges(:,3),GraphEdges(:,4),GraphEdges(:,5), ...
            'VariableNames',{'EndNodes' 'Width' 'Le' 'Lp'});
        NodeTable = table(GraphNodes(:,1),GraphNodes(:,2),GraphNodes(:,3),'VariableNames',{'X' 'Y','Type'});
        G = graph(EdgeTable,NodeTable);
        G = rmedge(G, 1:numnodes(G), 1:numnodes(G)); % Remove self loops
        G=rmnode(G,find(degree(G)==0));
    
    end    
end