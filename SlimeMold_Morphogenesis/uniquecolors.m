function [UC]=uniquecolors(ImagesFolder,imgType,mask)

srcFiles =dir([ImagesFolder '/' imgType]); % Gets info of the files inside the given folder with the given extension
nImages=length(srcFiles); %Number of images inside the folder 
UC=[];

for i =1:nImages %Loop through every image in the folder
    filename = strcat(ImagesFolder,'/',srcFiles(i).name); %Set the current Image identifier
    I = imread(filename); % Read the image 
    I_v=reshape(I,numel(I(:,:,1)),3); %Convert it to a vector (rows*columns,3)
    ListColors=I_v(mask,:); %Reduced vector (just non-mask positions) with 3 columns
    UC=unique([ListColors;UC],'rows');  
end

UC=unique(UC,'rows');
