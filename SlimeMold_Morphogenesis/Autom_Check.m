% Automated routine -  Calls all the set of routines for each one of the
% sub-folders under the root directory
clear 
close all
clc

% Get a list of all files and folders in this folder.
folder=uigetdir(pwd, 'Specify the folder containing the parent and result folders');
files = dir(folder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Print folder names to command window.

idcs   = strfind(folder,'\');
ResultsFolder=[folder(1:idcs(end)-1) filesep strcat('Results-',folder(idcs(end)+1:end))];
mkdir(ResultsFolder); %Create Results Folder

%LOOP to get masks
masks=cell(length(subFolders),1);
UC=cell(length(subFolders),1);
imgType = '*.tif'; % change based on image type
AutoMask=0;

for k = 3 : length(subFolders)
    masks{k}=getMask_auto(fullfile(folder,subFolders(k).name),imgType,AutoMask);
end

%LOOP to Classify and Analyze
for k = 3 : length(subFolders)
	fprintf('Sub folder #%d = %s\n', k, subFolders(k).name);
    ResultsPath=[ResultsFolder filesep strcat('Results-',subFolders(k).name)];
    mkdir(ResultsPath); %Create Results Folder

    tic
    
    ImagesFolder=fullfile(folder,subFolders(k).name);
    
    UC=uniquecolors(ImagesFolder,imgType,masks{k});

    funct_Trainer(ImagesFolder,ResultsPath,2,UC);
    
    funct_SMIA(ImagesFolder,masks{k},ResultsPath,subFolders(k).name);
    
    funct_PP(subFolders(k).name,ResultsPath)  
    
    close all
    toc
end
