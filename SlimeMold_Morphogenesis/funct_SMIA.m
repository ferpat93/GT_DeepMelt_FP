function []=funct_SMIA(ImagesFolder,mask,ResultsPath,Name)
warning off all
close all;

imgType = '*.tif'; % change based on image type
srcFiles =dir([ImagesFolder '/' imgType]); % Gets info of the files inside the given folder with the given extension
nImages=length(srcFiles); %Number of images inside the folder

%Get the classifier
Classifier = load(fullfile(ResultsPath,'Classifier'));
CM=Classifier.Classifier.ColorMap;

%FIRST IMAGE - SIZE AND MASK    
filename1 = fullfile(strcat(ImagesFolder,'/',srcFiles(1).name)); %get first image in the folder
I = imread(filename1); % open the image
I = I(:,:,1:3); %Erase transparency dimension
[rows, columns, ~] = size(I); %Get size of images
%mask=find(reshape(getMask(I),rows*columns,1));  % Get Mask to apply - vector of vector positions of non mask pixels
%npix=sum(sum(mask)); %num of image pixels inside region of interest

ImagesArray=zeros(rows*columns,nImages); % Cell array cotaining nImages number of mXnX3 matrices

%Classifier 
for i =1:nImages %Loop through every image in the folder
    filename = strcat(ImagesFolder,'/',srcFiles(i).name); %Set the current Image identifier
    I = imread(filename); % Read the image 
    Iv=reshape(I,numel(I(:,:,1)),3); %Convert it to a vector (rows*columns,3)
    I_nonMask=Iv(mask,:); %Reduced vector (just non-mask positions) with 3 columns
    
    if CM==2
        I_nonMask = rgb2lab(I_nonMask);
    elseif CM==3
        int=rgb2lab(I_nonMask);
        I_nonMask=int(:,2:3);
    end
    
    Identities=Classifier.Classifier.predictFcn(I_nonMask); %calls classifier, returns reduced binary vector 
    
    mappedI_v=4.*ones(length(Iv(:,1)),1); 
    mappedI_v(mask)=Identities; %Full vector (-1 for mask, 0 for non-SM, 1 for SM)
    
    ImagesArray(:,i)=mappedI_v;
end

Cfullpath = fullfile(ResultsPath,strcat(Name,'-ClassifiedImages')); % Location of ClassifiedImages
save(Cfullpath,'ImagesArray','rows','columns','-v7.3');



