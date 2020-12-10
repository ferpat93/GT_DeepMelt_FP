function [mask]=getMask_auto(ImagesFolder,imgType,AutoMask)

srcFiles =dir([ImagesFolder '/' imgType]); % Gets info of the files inside the given folder with the given extension

% Ask for mask 
filename1 = [ImagesFolder filesep srcFiles(1).name]; %get first image in the folder
I = imread(filename1); % open the image
I = I(:,:,1:3); %Erase transparency dimension
[rows, columns, ~] = size(I); %Get size of images
mask=find(reshape(getMask(I,AutoMask),rows*columns,1));  % Get Mask to apply - vector of vector positions of non mask pixels
