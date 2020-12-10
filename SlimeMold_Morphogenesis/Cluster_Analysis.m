function []=Cluster_Analysis(Name,Path)

struct = load(fullfile(Path,strcat(Name,'-ClassifiedImages')));

Images_BIN=struct.ImagesArray;
rows=struct.rows;
columns=struct.columns;
nImages=length(Images_BIN(1,:));

ShapeAnalysis=zeros(4,nImages);

%% Start Skeletonization 

%Images BW is a real binary image with SM=1 and 0 elsewhere
Images_BW=Images_BIN;
Images_BW(Images_BW==4)=0; % mask becomes void
Images_BW(Images_BW==2)=0; % redisuumm becomes void

rse=20; %Size of the dilation/erosion element
se = strel('disk', rse); % sets structural element as a disk of the given size

MinSize=0.005; % Min area of a region to be considered as a cluster - as percentage of total area 

for i=1:nImages
    % Remove interior pixels to leave an outline of the shapes.
    %BW2 = bwmorph(BW,'remove');
    
    Image=reshape(Images_BW(:,i),rows,columns); % Reshapes the binary image from column to matrix
    
    %figure(1)
    %imshow(Image)
    
    %   Perform filling operations
    
    % 'bridge' : Bridges unconnected pixels, that is, sets 0-valued pixels to 1 if they have two nonzero neighbors that are not connected. For example:
    Image=bwmorph(Image,'bridge'); 
    % 'fill' : Fills isolated interior pixels (individual 0s that are surrounded by 1s), such as the center pixel in this pattern.
    Image=bwmorph(Image,'fill');
    Image=imfill(Image,'holes'); % fills nonSM pixels inside a closed SM boundary
    % 'majority' : Sets a pixel to 1 if five or more pixels in its 3-by-3 neighborhood are 1s; otherwise, it sets the pixel to 0.
    Image=bwmorph(Image,'majority'); 
    % 'clean' : Removes isolated pixels (individual 1s that are surrounded by 0s), such as the center pixel in this pattern.
    Image=bwmorph(Image,'clean'); 
    
    %figure(2)
    %imshow(Image)
    %% Solidity and Circularity
    
    % Get biggest blob of SM e.g. biggest connected region
    %SingleBlob=zeros(rows,columns); % New matrix to store single largest blob
    %CCC=bwconncomp(Image); % Finds connected regions inside image
    %numPixels = cellfun(@numel,CCC.PixelIdxList); %tranform cell to mat array
    %[~,idx] = max(numPixels); % finds biggest connected region
    %SingleBlob(CCC.PixelIdxList{idx}) = 1; %stores the given connected region to te Single Blob array
    contour=bwconvhull(Image);
    %assignin('base','contour',contour)
    %figure(3)
    %imshow(contour)
    %ShapeIndexes = regionprops(SingleBlob,'Area','Perimeter','Eccentricity','Solidity'); % Gets its properties in a struct
    ShapeIndexes = regionprops(contour,'Area','Perimeter','Eccentricity'); % Gets its properties in a struct
    ShapeIndexes.Solidity = ((sum(sum(Image)))/(sum(sum(contour)))); % Ratio between real area and convex area
    ShapeIndexes.Circularity = (ShapeIndexes.Perimeter .^ 2) ./ (4 * pi * ShapeIndexes.Area); % Circularity := (Perimeter .^ 2) ./ (4 * pi * area);
    %assignin('base','ShapeIndexes',ShapeIndexes)
    %% Number of Clusters
    Eroded_I = imopen(Image, se); % Use structure to erode and disconnect entities
    CC=bwconncomp(Eroded_I);
    %assignin('base','CC_eroded',CC)
    %figure(4)
    %imshow(Eroded_I)
    
    numPixels = cellfun(@numel,CC.PixelIdxList); %tranform cell to mat array
    SizeIndex = numPixels>MinSize*rows*columns; % finds biggest connected region
    ShapeIndexes.NumClusters=sum(SizeIndex);
    ShapeAnalysis(:,i)=[ShapeIndexes.Eccentricity;ShapeIndexes.Solidity;ShapeIndexes.Circularity;ShapeIndexes.NumClusters];
    %assignin('base','ShapeIndexes',ShapeIndexes)
end

Cfullpath = fullfile(Path,strcat(Name,'-ShapeAnalysis')); % Gets info of the files inside the given folder with the given extension
save(Cfullpath,'ShapeAnalysis','-v7.3');
