% Function Get Network
% Input Image - Output Graph

function [labeled_clusters,Veins,G,centroids] = getNetwork_Fusion(ind_im,Io,mask,centroids)

%% Preprocess image -> Raw to Full Binary

    % pixels per mm
    D_dish = 90; %mm
    ppmm = size(Io,1)/D_dish;
    
    %%% Enhance contrast - grayscale - binarize

    % Enhance
    sigma = 0.4;
    alpha = 0.675;
    I = locallapfilt(Io, sigma, alpha); % Laplacian filtering -  
    %I = locallapfilt(Io, sigma, alpha, 'NumIntensityLevels', 20); % Laplacian filtering -  

    lab = rgb2lab(I); % to LAB
    N=lab(:,:,1)-lab(:,:,3); % Compound 
    GI=(N-min(N,[],'all'))/(max(N,[],'all')-min(N,[],'all')); % Normalize
    %GI = lab(:,:,1)./100; % rescale L component (equivalent to greyscale)
    %GI = rgb2gray(I)./255;
    %GI = adapthisteq(GI); % Optional contrast enhancement
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
            pause
        else
            centroids=numPixels.Centroid(idx,:);
        end
    end 
    
    nCells = size(centroids,1); % 1 or 2 initial cells
    
    %% Split cells in two
    
    % Connected Print %
    BW = ConnectedPrint(BW,T);   
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
    
    %%% Identify Clusters and Measure Vein Width
    
    maxVeinSize = round(1*ppmm)+1; % 4mm meaning 8mm min dimension.
    labeled_clusters=zeros(size(BW));
    
    for cell=1:NC
         
        seE = strel('disk', double(maxVeinSize)); % erode and dilate image
        EI = imerode(BW(:,:,cell),seE);

        seD = strel('disk', (round(1.5*maxVeinSize)+1));
        DI = imdilate(EI,seD);

        clustersI = bwlabel(DI.*(BW(:,:,cell))); % Get Clusters

        clusters_label = unique(clustersI.*EI);
        clustersI = imfill(ismember(clustersI, clusters_label(2:end)),'holes');
        labeled_clusters(:,:,cell) = int8(bwlabel(clustersI));
        %nClusters = max(labeled_clusters(:));

        Full_BI = (BW(:,:,cell) + clustersI)>0; % after filling inside of cluster
        Idist = bwdist(~Full_BI); % Binary distance (without clusters!)
        Veins(:,:,cell) = Idist.*(~clustersI);
        
        %figure()
        %imagesc(Idist)

        % ---- Uncomment if graph is needed --- %
        skel_refined = bwskel(Full_BI>0,'MinBranchLength',10);
        skel_dist = skel_refined.*Idist;
        G{cell} = BinToGraph(skel_dist,labeled_clusters(:,:,cell));        
        %G = [];

        % OJO
        %figure()
        %imshow(skel_refined)
        %title('Skeleton')
      
    end
    
    %labeled_clusters = labeled_clusters(:,:,1)-labeled_clusters(:,:,2);
    %Veins = Veins(:,:,1)-Veins(:,:,2);
    
    %% NOT ACTIVE! : lines that deal with more than 1 CC -> Not accurate! Can be disconnected by cluster
    %CC = bwconncomp(veins);
    %numPixels = cellfun(@numel,CC.PixelIdxList);
    %[~,idx] = max(numPixels);
    %veins(:) = 0;
    %veins(CC.PixelIdxList{idx}) = 1;


end

%% Auxiliar Functions
