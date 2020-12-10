% Check Progress of analyses
clear 
close all
clc

Output={'Classifier.mat','-ClassifiedImages.mat','-Evolution_over_time.mat','-Spatial_analysis.mat','-BiDish.mat','-ShapeAnalysis.mat'};
Functions={'funct_Trainer(ImagesFolder,ResultsPath,2,uniquecolors(ImagesFolder,imgType,masks{k}))',...
    'funct_SMIA(ImagesFolder,masks{k},ResultsPath,Name)',...
    'funct_PP(Name,Specific_ResultsPath)',...
    'FunctSpatialAnalysis(circles,PhotosFolder,Name,Specific_ResultsPath,k)',...
    'BiDish(circles{k-2,1},Specific_ResultsPath,Name)',...
    'Cluster_Analysis(Name,Specific_ResultsPath)'};

load('C:\Users\lfp3\OneDrive - Georgia Institute of Technology\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\circles.mat')
imgType = '*.tif'; % change based on image type

% Get a list of all files and folders in this folder.
BaseFolder=uigetdir(pwd, 'Specify the folder containing the PHOTOS and RESULTS folders for the analysis');
PhotosFolder=fullfile(BaseFolder,'PHOTOS');
ResultsFolder=fullfile(BaseFolder,'RESULTS');

files = dir(PhotosFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

Ledger=zeros(length(subFolders)-2,5);

%LOOP to get masks
%masks=cell(length(subFolders),1);
%UC=cell(length(subFolders),1);
%AutoMask=0;

%for k = 3 : length(subFolders)
%    masks{k}=getMask_auto(fullfile(folder,subFolders(k).name),imgType,AutoMask);
%end


%LOOP to Classify and Analyze
for k = 3 : length(subFolders)
   
    fprintf('Sub folder #%d = %s\n', k-2, subFolders(k).name);
    Name=subFolders(k).name;
    Specific_ResultsPath=[ResultsFolder filesep 'Results-' subFolders(k).name];
    ImagesFolder=[PhotosFolder filesep subFolders(k).name];
    
    for f=2:numel(Output)
        nameOut=fullfile(Specific_ResultsPath,strcat(Name,Output{f}));
        close all
        if (exist(nameOut,'file'))>0 
            Ledger(k-2,f)=1;
        else
            funct=Functions{f};
            eval(funct)
        end
    
    end
    
    close all
 
end


