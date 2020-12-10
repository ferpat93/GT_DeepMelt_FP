function [UC,mask]=uniquecolors_wmask(ImagesFolder,imgType,AutoMask)

srcFiles =dir([ImagesFolder '/' imgType]); % Gets info of the files inside the given folder with the given extension
nImages=length(srcFiles); %Number of images inside the folder 
UC=[];

% Ask for mask 
filename1 = fullfile(strcat(ImagesFolder,'/',srcFiles(1).name)); %get first image in the folder
I = imread(filename1); % open the image
I = I(:,:,1:3); %Erase transparency dimension
[rows, columns, ~] = size(I); %Get size of images
mask=find(reshape(getMask(I,AutoMask),rows*columns,1));  % Get Mask to apply - vector of vector positions of non mask pixels

for i =1:nImages %Loop through every image in the folder
    filename = strcat(ImagesFolder,'/',srcFiles(i).name); %Set the current Image identifier
    I = imread(filename); % Read the image 
    I_v=reshape(I,numel(I(:,:,1)),3); %Convert it to a vector (rows*columns,3)
    ListColors=I_v(mask,:); %Reduced vector (just non-mask positions) with 3 columns
    UC=unique([ListColors;UC],'rows');  
end

UC=unique(UC,'rows');
