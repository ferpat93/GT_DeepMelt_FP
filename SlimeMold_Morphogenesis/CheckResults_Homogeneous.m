% Check Progress of analyses
clear 
close all
clc

Output={'Classifier.mat','-ClassifiedImages.mat','-Evolution_over_time.mat','-ShapeAnalysis.mat'};
Functions={'funct_Trainer([ImagesFolder filesep srcFiles(1).name],Specific_ResultsPath,3,uniquecolors([ImagesFolder filesep srcFiles(1).name],imgType,masks{k}))',...
    'funct_SMIA([ImagesFolder filesep srcFiles(1).name],masks{k},Specific_ResultsPath,Name)',...
    'funct_PP(Name,Specific_ResultsPath)',...
    'Cluster_Analysis(Name,Specific_ResultsPath)'};

load('C:\Users\lfp3\OneDrive - Georgia Institute of Technology\Slime-Mold-GT-Tolouse\SPOT-ANALYSIS\RESULTS\circles.mat')
imgType = '*.tif'; % change based on image type

% Get a list of all files and folders in this folder.
BaseFolder=uigetdir('C:\Users\lfp3\OneDrive - Georgia Institute of Technology\Slime-Mold-GT-Tolouse\HOMOGENEOUS-ANALYSIS', 'Specify the folder containing the PHOTOS and RESULTS folders for the analysis');
PhotosFolder=fullfile(BaseFolder,'PHOTOS');
ResultsFolder=fullfile(BaseFolder,'RESULTS');

files = dir(PhotosFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

Ledger=zeros(length(subFolders)-2,5);

%masksPath = fullfile(ResultsFolder,'Masks'); % fullpath of the masks
   
%if (exist(masksPath,'file'))>0 
%    load(masksPath);
%else
    %%LOOP to get masks
    %%for k = 3 : length(subFolders)
     %   masks{k}=getMask_auto(fullfile(PhotosFolder,subFolders(k).name),imgType,0);
    %end
    %save(masksPath,'masks','-v7.3');
%end

%LOOP to Classify and Analyze
for k = 29 : length(subFolders)
   
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


