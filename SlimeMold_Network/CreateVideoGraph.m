% Create Video - Graph

close all
% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.065 0.065], [0.2 0.05], [0.075 0.04]); % subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
% if ~make_it_tight,  clear subplot;  end

DataFolder = 'D:\Slime_Mold_Network\Data\Fusion';

Folder = '24_12_19_Chamber_2';   
%Folder = '24_12_19_Chamber_1';

ROI_filename = 'ROI_file.csv';
srcFiles_filename = 'srcFiles.mat';
ImgType = '*.jpg';

dishes = [1];
t_int = 20;

ROI = csvread(fullfile(DataFolder,Folder,ROI_filename));      
nE = size(ROI,1);
    
srcFiles = load(fullfile(DataFolder,Folder,srcFiles_filename));
srcFiles = extractfield(srcFiles,'srcFiles'); srcFiles = srcFiles{1};

times = 1:t_int:size(srcFiles,1);
% Video Parameters
%
filepathE = fullfile(pwd,'Video_Entities');

FPS=2; % Frames per second of video
nF=numel(times); %Total number of frames to display
F(1) = struct('cdata',[],'colormap',[]); % Creates struct variable to store figure frames
%

for i=1:numel(dishes)
    e = dishes(i);
    mask = uint8(createCircleMask(ROI(e,3)));    
    centroids = [ 10 500; 300 500];
    
    
    for ti=1:numel(times)
        time = times(ti);
        path = fullfile(DataFolder,Folder,srcFiles.name{time});

        [BW,I]= getBinary(path,ROI(e,:),centroids);
        [labeled_clusters,skeleton]=GetClusters_Skeleton(BW);

        B = mask.*superposedImage(I,BW,skeleton,labeled_clusters);

        figure(1)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        axis equal % Sets equal scale for axis
       
        imshow(B)
        drawnow %Update plot
    
        F(ti) = getframe(gcf); %Get frame and store it  
    
    end
end

v = VideoWriter(filepathE); % Write video in the desired location
v.Quality=100;
v.FrameRate = round(FPS); %Sets frames per second
open(v) %Activate video
writeVideo(v,F) %Save it
close(v) %Deactivate video

%% Auxiliar Functions

function [B] = superposedImage(I,SM,skeleton,clusters_labeled)
    
    %border = bwperim(SM);
    clusters = clusters_labeled>0;
    
    colors = [0.9 0.67 0; 0 0.4 0.4; 0.5 0.0 0];
    position = [1 1];

    A = cat(3,3.*skeleton,2.*clusters,SM);
    L = max(A,[],3);
    B = labeloverlay(I,L,'Colormap',colors);
    %B = uint8(insertText(B,position,'Control - An = 4','BoxOpacity',1,'BoxColor','w'));
    
end

function [labeled_clusters,skel_refined] = GetClusters_Skeleton(BW) % returns a table

    % Split BW into Veins/Clusters
        
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
    labeled_clusters = bwlabel(clustersI);
    Full_BI = (BW + clustersI)>0; % after filling inside of clusters
    skel_refined = bwskel(Full_BI>0,'MinBranchLength',10);

end

function [mask] = createCircleMask(L)
    
    [xx,yy] = meshgrid(1:L,1:L);
    mask = false(L,L);
    mask = mask | hypot(xx - L/2, yy - L/2) <= 0.99*(L/2);

end

function [BW,I] = getBinary(path,ROI,centroids)

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

