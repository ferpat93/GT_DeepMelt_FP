% Automated routine -  Calls all the set of routines for each one of the
% sub-folders under the root directory

close all
clc

% Get a list of all files and folders in this folder.
folder=uigetdir(pwd, 'Specify the folder containing the Image folders for the analysis');
files = dir(folder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.

idcs   = strfind(folder,'\');
ResultsFolder=[folder(1:idcs(end)-1) filesep 'RESULTS']; % ojo something weird
mkdir(ResultsFolder); %Create Results Folder

%LOOP to get masks

UC=cell(length(subFolders),1);
imgType = '*.tif'; % change based on image type
AutoMask=0; % 0 for requesting new, % 1 for attempting automatic, % 2 for saved

%{
masks=cell(length(subFolders),1);
for k = 3 : length(subFolders)
    masks{k}=getMask_auto(fullfile(folder,subFolders(k).name),imgType,AutoMask);
end
%}


%LOOP to Classify and Analyze
for k = 10 : 10%length(subFolders)
	fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
    ResultsPath=[ResultsFolder filesep strcat('Results-',subFolders(k).name)];
    mkdir(ResultsPath); %Create Results Folder

    tic
    
    ImagesFolder=fullfile(folder,subFolders(k).name);
    
    %UC=uniquecolors(ImagesFolder,imgType,masks{k});

    %funct_Trainer(ImagesFolder,ResultsPath,3,UC);
    %funct_Trainer_k(ImagesFolder,ResultsPath,3,UC);
    
    %funct_SMIA(ImagesFolder,masks{k},ResultsPath,subFolders(k).name);
    
    %funct_PP(subFolders(k).name,ResultsPath) 
    
    imgDiff=funct_HDKM(ImagesFolder,ResultsPath,subFolders(k).name);
    
    %Cluster_Analysis(subFolders(k).name,ResultsPath)
    
    close all
    toc
end
