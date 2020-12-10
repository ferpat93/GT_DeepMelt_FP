function [Images_DIF]=funct_HDKM(ImagesPath,ResultsPath,Name)
% High Dimensional K-means
warning off all
close all;

nV=1; %number of pixels to the side of window center;

imgType = '*.tif'; % change based on image type
srcFiles =dir([ImagesPath '/' imgType]); % Gets info of the files inside the given folder with the given extension

struct = load(fullfile(ResultsPath,strcat(Name,'-ClassifiedImages')));
Images_TRIN=struct.ImagesArray;
rows=struct.rows;
columns=struct.columns;
nImages=length(Images_TRIN(1,:));

%FIRST IMAGE - SIZE AND MASK    
filename1 = fullfile(strcat(ImagesPath,'/',srcFiles(1).name)); %get first image in the folder
Io = imread(filename1); % open the image
Io = reshape(sum(Io(:,:,1:3),3),rows*columns,1); %Erase transparency dimension

Images_DIF=Images_TRIN;
Images_NH=nan(rows*columns,nImages);
Data=[];
for i =1:nImages %Loop through every image in the folder
    
    %BackGround Substraction
    
    filename = strcat(ImagesPath,'/',srcFiles(i).name); %Set the current Image identifier
    I = imread(filename); % Read the image 
    Iv=reshape(sum(I,3),numel(I(:,:,1)),1); %Convert it to a vector (rows*columns,3)
    Images_DIF(:,i)=abs(Iv-Io);
    
    %Neighborhood Analysis:    
    %lab=reshape(rgb2lab(I),numel(I(:,:,1)),3);
    
    %Data=[Data;[Images_DIF(:,i), lab(:,2:3)]];
     
    %Images_NH(:,i)=NeighborCount(reshape(Images_TRIN(:,i),rows,columns),nV);
    
end

end


function [NC]=NeighborCount(I,nN)

rows=size(I,1);
columns=size(I,2);


% Value Codes
% SM = 1
% Agar = 2
% Blob = 3
% Mask = 4

fun = @(x) median(x(:));

NC = nlfilter(I,[1+2*nN 1+2*nN],fun); 
NC(NC==4) = I(NC==4);

NC=reshape(NC,columns*rows,1);

end

